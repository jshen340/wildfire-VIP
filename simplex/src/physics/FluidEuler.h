//////////////////////////////////////////////////////////////////////////
// Incompressible fluid solver on an Eulerian grid
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "MacGrid.h"
#include "Advection.h"
#include "BoundaryCondition.h"
#include "ProjectionMeta.h"
#include "ProjectionCPX.h"
#include "KrylovSolver.h"
#include "Hashtable.h"
#include "Timer.h"
#include "Interpolation.h"
#include "FluidFunc.h"

template<int d> class FluidEuler
{Typedef_VectorDii(d);
public:
	//data
	MacGrid<d> mac_grid;
	FaceField<real,d> velocity;
	FaceField<real,d> alpha;
	Field<ushort,d> type;		////supporting type: Fluid, Solid
	BoundaryConditionMacGrid<d> bc=BoundaryConditionMacGrid<d>(mac_grid);
	std::conditional_t<d == 2, Field<real, 2>, Field<VectorD, 3> > vorticity;
	//projection
	std::shared_ptr<ProjectionCPX<real, d>> projection;
	//settings
protected:
	bool use_body_force=false;
	VectorD g=VectorD::Unit(1)*(real)-9.8;
	bool use_maccormack = true;
	bool use_zero_extrapolation = true; //zero extrapolation for pos out of field during advection  
public:
	FluidEuler(json& j) {
		projection = std::make_shared<ProjectionCPX<real, d>>(
			j.at("projection"),
			&mac_grid, &type, &bc
			);
		use_body_force = j.value<bool>("use_body_force", false);
	}

	virtual void Initialize(const VectorDi& cell_counts,const real dx,const VectorD& domain_min=VectorD::Zero())
	{
		mac_grid.Initialize(cell_counts,dx,domain_min);
		Initialize_Fields();
	}

	virtual void Initialize_Fields()
	{
		velocity.Resize(mac_grid.grid.cell_counts,(real)0);
		alpha.Resize(mac_grid.grid.cell_counts,(real)1);
		type.Resize(mac_grid.grid.cell_counts,(ushort)CellType::Fluid);
		vorticity.Resize(mac_grid.grid.cell_counts);
	}

	virtual void Advance(const real dt)
	{
		Timer timer;timer.Reset();
		Advection(dt);
		timer.Elapse_And_Output_And_Reset("Advection");
		Apply_Body_Forces(dt);
		Enforce_Incompressibility(dt);
		timer.Elapse_And_Output_And_Reset("Projection");
	}

	virtual void Advection(const real dt)
	{
		MacGrid<d> ghost_grid=mac_grid;FaceField<real,d> ghost_velocity=velocity;
		if(use_maccormack){ Advection::MacCormack(dt, ghost_grid, ghost_velocity, velocity, use_zero_extrapolation,true); }
		else{ Advection::Semi_Lagrangian(dt, ghost_velocity, ghost_grid, ghost_velocity, mac_grid, velocity, use_zero_extrapolation); }
		Update_Cell_Types();
	}

	virtual void Update_Cell_Types()
	{
		int cell_num=mac_grid.grid.cell_counts.prod();
		#pragma omp parallel for
		for(int i=0;i<cell_num;i++){VectorDi cell=mac_grid.grid.Cell_Coord(i);
			type(cell)=bc.Is_Psi_D(cell)?(ushort)CellType::Solid:(ushort)CellType::Fluid;}
	}

	virtual void Apply_Body_Forces(const real dt)
	{
		if (!use_body_force){return;}
		for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
			#pragma omp parallel for
			for(int i=0;i<face_num;i++){
				VectorDi face=mac_grid.Face_Coord(axis,i);
				velocity.face_fields[axis](face)+=g[axis]*dt;}}
	}

	virtual void Enforce_Incompressibility(const real dt)
	{
		Enforce_Boundary_Conditions();
		projection->Project(dt, velocity);
	}

	virtual void Enforce_Boundary_Conditions()
	{
		for(auto p:bc.psi_N_values){int axis=p.first[0];int face_index=p.first[1];real value=p.second;
			velocity.face_fields[axis].array[face_index]=value;}
	}

	//////////////////////////////////////////////////////////////////////////
	////helper functions

	real Max_Abs(const FaceField<real,d>& field_q) const
	{
		real max_abs=(real)0;
		for (int i = 0;i < d;i++) {
			int n = (int)field_q.face_fields[i].array.size();
			for (int j = 0;j < n;j++) {
				real abs_v = abs(field_q.face_fields[i].array[j]);
				if (abs_v > max_abs)max_abs = abs_v;
			}
		}
		return max_abs;
	}

	real Kinetic_Energy(const FaceField<real, d>& field_q) const 
	{
		real kinetic_energy = 0;

		int cell_num = field_q.mac_grid.grid.cell_counts.prod();
		#pragma omp parallel for reduction(+:kinetic_energy)
		for (int i = 0; i < cell_num; i++) {
			const VectorDi cell = mac_grid.grid.Cell_Coord(i);
			VectorD cell_v = VectorD::Zero();
			for (int axis = 0; axis < d; axis++) {
				VectorDi left_face = mac_grid.Cell_Left_Face(axis, cell);
				VectorDi right_face = mac_grid.Cell_Right_Face(axis, cell);
				cell_v[axis] = (real)0.5 * (field_q(axis, left_face) + field_q(axis, right_face));
			}
			kinetic_energy += cell_v.dot(cell_v);
		}
		return (real)0.5 * kinetic_energy * pow(mac_grid.grid.dx, d);
	}

	real Enstrophy()
	{
		FluidFunc::Curl_On_Cell(mac_grid, velocity, vorticity);
		real enstrophy = 0;
		int cell_num = mac_grid.grid.cell_counts.prod();
		#pragma omp parallel for reduction(+:enstrophy)
		for (int i = 0; i < cell_num; i++) {
			const VectorDi cell = mac_grid.grid.Cell_Coord(i);
			if constexpr (d == 2) { enstrophy += pow(vorticity(cell), 2); }
			else { enstrophy += vorticity(cell).dot(vorticity(cell)); }
		}

		return (real)0.5 * enstrophy * pow(mac_grid.grid.dx, d);
	}

	void Divergence() {
		real total_divergence = 0;
		for (int i = 0; i < mac_grid.grid.cell_counts.prod(); i++) {
			VectorDi cell = mac_grid.grid.Cell_Coord(i);
			real divergence = 0;
			for (int axis = 0; axis < d; axis++) {
				divergence += velocity(axis, mac_grid.Cell_Right_Face(axis, cell));
				divergence -= velocity(axis,mac_grid.Cell_Left_Face(axis, cell));
			}
			divergence /= mac_grid.grid.dx;
			std::cout << "Divergence for cell " << cell.transpose() << " is: " << divergence << std::endl;
			total_divergence += divergence;
		}
		std::cout << "Total divergenec is: " << total_divergence << std::endl;
	}
};