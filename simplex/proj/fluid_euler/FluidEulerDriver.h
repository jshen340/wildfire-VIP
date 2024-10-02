//////////////////////////////////////////////////////////////////////////
// Euler fluid driver
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __FluidEulerDriver_h__
#define __FluidEulerDriver_h__
#include "Common.h"
#include "File.h"
#include "GeometryPrimitives.h"
#include "Particles.h"
#include "MarchingCubes.h"
#include "Driver.h"
#include "Poisson.h"
#include "FluidEuler.h"
#include "RenderFunc.h"

template<int d> class FluidEulerDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;
public:
	std::shared_ptr<FluidEuler<d>> fluid;

	FluidEulerDriver() {}
	//FluidEulerDriver():fluid(SolverType::CPX_GPU){}

	real CFL() const
	{
		//real epsilon=(real)1e-5;
		//return cfl*fluid.mac_grid.grid.dx/(fluid.Max_Abs(fluid.velocity)+epsilon);
		return (real)0.025;
	}

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		fluid->Advance(dt);
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);
		const auto mac_grid = fluid->mac_grid;
		if(frame==0){
			std::string file_name=frame_dir+"/grid";
			mac_grid.grid.Write_To_File_3d(file_name);
			std::cout<<"Write to file "<<file_name<<std::endl;}
		
		////Write velocity
		fluid->velocity.Write_To_File_3d(frame_dir + "/velocity");
		////Write BC
		RenderFunc::Write_Dirichlet_Boundary_Conditions(frame_dir + "/psi_D", mac_grid, fluid->bc);
		RenderFunc::Write_Neumann_Boundary_Conditions(frame_dir + "/psi_N", mac_grid, fluid->bc);

		////Write fluid type
		{std::string file_name=frame_dir+"/fluid_type";
		Field<real,d> fluid_type;fluid_type.Resize(mac_grid.grid.cell_counts);
		iterate_cell(iter, mac_grid.grid){const VectorDi& cell=iter.Coord();
			if(fluid->type(cell)==(ushort)CellType::Fluid)fluid_type(cell)=(real)0;
			else fluid_type(cell)=(real)1;}
		fluid_type.Write_To_File_3d(file_name);}

		////Write energy
		if constexpr(d==2){
		real energy = fluid->Kinetic_Energy(fluid->velocity);
		real enstrophy;
		enstrophy = fluid->Enstrophy();
		std::cout << "Kinetic Energy for frame " << frame << ": " << energy << std::endl;
		std::cout << "Enstrophy for frame " << frame << ": " << enstrophy << std::endl;
		std::string file_name = output_dir + "/0/energy.txt";
		std::string info = std::to_string(frame) + "\t" + std::to_string(energy) + "\t" + std::to_string(enstrophy) + "\n";
		if (frame == 0) { File::Write_Text_To_File(file_name, "Frame\tenergy\tenstrophy\n"); }
		File::Append_Text_To_File(file_name, info); }

		////Write vorticity
		if constexpr (d == 2) {
			std::string file_name = frame_dir + "/vorticity";
			fluid->vorticity.Write_To_File_3d(file_name);
		}

		std::cout<<"Write to frame "<<frame<<std::endl;
	}

	virtual void Update_Face_Fraction(LevelSet<d>& _ls_solid)
	{
		const auto mac_grid = fluid->mac_grid;
		iterate_face(axis, iter, mac_grid) {
			const VectorDi& face = iter.Coord();
			VectorDi node_0 = mac_grid.Face_Incident_Node(axis, face, 0);
			VectorDi node_1 = mac_grid.Face_Incident_Node(axis, face, 1);
			real phi_0 = _ls_solid.Phi(mac_grid.grid.Node(node_0));
			real phi_1 = _ls_solid.Phi(mac_grid.grid.Node(node_1));
			// no need to multiply dx
			if (phi_0 < (real)0 && phi_1 >(real)0) {
				fluid->alpha.face_fields[axis](face) = phi_0 / (phi_0 - phi_1);
			}
			else if (phi_0 > (real)0 && phi_1 < (real)0) {
				fluid->alpha.face_fields[axis](face) = phi_1 / (phi_1 - phi_0);
			}
			else if (phi_0 > (real)0 && phi_1 > (real)0) {
				fluid->alpha.face_fields[axis](face) = (real)0;
			}
			else {
				fluid->alpha.face_fields[axis](face) = (real)1;
			}
		}
	}
	
	virtual void Initialize()
	{
		int s = scale; real length = (real)1; VectorDi cell_counts = VectorDi::Ones(); cell_counts[0] = s;
		cell_counts[1] = cell_counts[0] / 2; if (d > 2)cell_counts[2] = cell_counts[0] / 2;

		switch (test) {
		case 1: {	////Karman vortex street
			json j_projection; j_projection = { {"update_A",false} };
			json j_fluid = { { "projection", j_projection } };
			fluid = std::make_shared<FluidEuler<d>>(j_fluid);
			fluid->Initialize(cell_counts, (real)length / cell_counts[0]);
			fluid->velocity.Fill((real).2, 0);

			const MacGrid<d> mac_grid = fluid->mac_grid;
			////Source
			real source_speed = (real)2; int axis = 0;
			iterate_face_in_one_dim(axis, iter, mac_grid) {
				const VectorDi& face = iter.Coord();
				if (face[0] == 0/*||face[0]==fluid.mac_grid.face_grids[0].node_counts[0]-1*/)fluid->bc.Set_Psi_N(0, face, source_speed);
			}
			////Wall
			{axis = 1; iterate_face_in_one_dim(axis, iter, mac_grid) {
				const VectorDi& face = iter.Coord();
				if (face[axis] == 0 || face[axis] == mac_grid.face_grids[axis].node_counts[axis] - 1)fluid->bc.Set_Psi_N(axis, face, 0);
			}}
			if (d == 3) {
				axis = 2; iterate_face_in_one_dim(axis, iter, mac_grid) {
					const VectorDi& face = iter.Coord();
					if (face[axis] == 0 || face[axis] == mac_grid.face_grids[axis].node_counts[axis] - 1)fluid->bc.Set_Psi_N(axis, face, 0);
				}
			}
			////Solid
			VectorD box_center = (real).5 * (mac_grid.grid.domain_min + mac_grid.grid.domain_max); box_center[0] = (real).2;
			//VectorD box_size=VectorD::Ones()*(real).025;Box<d> obstacle(box_center-box_size,box_center+box_size);//
			real r = (real).04; Sphere<d> obstacle(box_center, r);
			iterate_cell(iter, mac_grid.grid) {
				const VectorDi& cell = iter.Coord(); const VectorD& pos = mac_grid.grid.Center(cell);
				if (obstacle.Inside(pos)) {
					fluid->bc.Set_Psi_D(cell, (ushort)CellType::Solid);
					for (int axis = 0; axis < d; axis++)for (int side = 0; side < 2; side++) { VectorDi face = mac_grid.Cell_Incident_Face(axis, cell, side); fluid->bc.Set_Psi_N(axis, face, (real)0); }
				}
			}
		}break;
		case 2: {	////Karman vortex street with accurate solid boundary
			json j_projection = json(json::value_t::object);//empty object to prevent parsing error
			json j_fluid = { { "projection", j_projection } };
			fluid = std::make_shared<FluidEuler<d>>(j_fluid);
			fluid->Initialize(cell_counts, (real)length / cell_counts[0]);
			fluid->velocity.Fill((real).2, 0);

			const auto mac_grid = fluid->mac_grid;
			////Source
			real source_speed = (real)2; int axis = 0;
			iterate_face_in_one_dim(axis, iter, mac_grid) {
				const VectorDi& face = iter.Coord();
				if (face[0] == 0/*||face[0]==fluid.mac_grid.face_grids[0].node_counts[0]-1*/)fluid->bc.Set_Psi_N(0, face, source_speed);}
			////Wall
			axis = 1; iterate_face_in_one_dim(axis, iter, mac_grid) {
				const VectorDi& face = iter.Coord();
				if (face[axis] == 0 || face[axis] == mac_grid.face_grids[axis].node_counts[axis] - 1)fluid->bc.Set_Psi_N(axis, face, 0);}
			////Solid
			VectorD box_center = (real).5 * (mac_grid.grid.domain_min + mac_grid.grid.domain_max); box_center[0] = (real).2;
			//VectorD box_size=VectorD::Ones()*(real).025;Box<d> obstacle(box_center-box_size,box_center+box_size);//
			real r = (real).04; Sphere<d> obstacle(box_center, r);
			iterate_cell(iter, mac_grid.grid) {
				const VectorDi& cell = iter.Coord(); const VectorD& pos = mac_grid.grid.Center(cell);
				if (obstacle.Inside(pos)) {
					fluid->bc.Set_Psi_D(cell, (ushort)CellType::Solid);
					for (int axis = 0; axis < d; axis++)for (int side = 0; side < 2; side++) { VectorDi face = mac_grid.Cell_Incident_Face(axis, cell, side); fluid->bc.Set_Psi_N(axis, face, (real)0); }}}

			// solid boundary levelset
			LevelSet<d> ls_solid; ls_solid.Initialize(mac_grid.grid);
			iterate_cell(iter, mac_grid.grid) {
				const VectorDi& cell = iter.Coord();
				VectorD pos = mac_grid.grid.Center(cell);
				ls_solid.phi(cell) = -obstacle.Phi(pos);}
			// update Alpha from solid levelset
			Update_Face_Fraction(ls_solid);
		}break;
		//case 3:{	////smoke plume
		//	fluid.Initialize(cell_counts,(real)length/cell_counts[0]);
		//	fluid.velocity.Fill((real).1,0);

		//	//real source_speed=(real)1;
		//	VectorD box_center=fluid.mac_grid.grid.Center();

		//	box_center[0]=(real).2;VectorD box_size=VectorD::Ones()*(real).025;box_size[0]=fluid.mac_grid.grid.dx*(real)2;
		//	Box<d> box(box_center-box_size,box_center+box_size);

		//	////Source
		//	//int axis=0;
		//	//iterate_face_in_one_dim(axis,iter,fluid.mac_grid){const VectorDi& face=iter.Coord();
		//	//	if(face[0]==0/*||face[0]==fluid.mac_grid.face_grids[0].node_counts[0]-1*/)fluid.bc.Set_Psi_N(0,face,source_speed);}

		//	////Initialize source
		//	VectorD source_velocity=AuxFunc::V<d>((real)1.,(real).1,(real).2);
		//	iterate_cell(iter,fluid.mac_grid.grid){const VectorDi& cell=iter.Coord();
		//		const VectorD& pos=fluid.mac_grid.grid.Center(cell);
		//			if(box.Inside(pos)){fluid.bc.Set_Psi_D(cell,(ushort)CellType::Solid);
		//				for(int axis=0;axis<d;axis++)
		//					for(int side=0;side<2;side++){
		//					VectorDi face=fluid.mac_grid.Cell_Incident_Face(axis,cell,side);fluid.bc.Set_Psi_N(axis,face,source_velocity[axis]);}}}

		//	////Solid
		//	VectorD solid_center=(real).5*(fluid.mac_grid.grid.domain_min+fluid.mac_grid.grid.domain_max);solid_center[0]=(real).4;
		//	real r=(real).04;Sphere<d> obstacle(solid_center,r);
		//	iterate_cell(iter,fluid.mac_grid.grid){const VectorDi& cell=iter.Coord();const VectorD& pos=fluid.mac_grid.grid.Center(cell);
		//		if(obstacle.Inside(pos)){fluid.bc.Set_Psi_D(cell,(ushort)CellType::Solid);
		//			for(int axis=0;axis<d;axis++)for(int side=0;side<2;side++){VectorDi face=fluid.mac_grid.Cell_Incident_Face(axis,cell,side);fluid.bc.Set_Psi_N(axis,face,(real)0);}}}

		//	fluid.projection.update_A=false;
		//}break;
		//case 4:{	////lid-driven cavity
		//	cell_counts=VectorDi::Ones()*s;
		//	fluid.Initialize(cell_counts,(real)length/cell_counts[0]);
		//	fluid.velocity.Fill((real).2,0);

		//	real source_speed=(real)2;
		//	int axis=0;
		//	iterate_face_in_one_dim(axis,iter,fluid.mac_grid){const VectorDi& face=iter.Coord();
		//		if(face[1]==0)fluid.bc.Set_Psi_N(0,face,0);	////bottom
		//		else if(face[1]==fluid.mac_grid.face_grids[axis].node_counts[1]-1)fluid.bc.Set_Psi_N(0,face,source_speed);	////top
		//		else if(face[0]<1||face[0]>fluid.mac_grid.face_grids[axis].node_counts[0]-2)fluid.bc.Set_Psi_N(axis,face,0);}	////left and right

		//	axis=1;
		//	iterate_face_in_one_dim(axis,iter,fluid.mac_grid){const VectorDi& face=iter.Coord();
		//		if(face[0]==0||face[0]==fluid.mac_grid.face_grids[axis].node_counts[0]-1||	////left and right
		//			face[1]<1||face[1]>fluid.mac_grid.face_grids[axis].node_counts[1]-2)fluid.bc.Set_Psi_N(axis,face,0);}

		//	iterate_cell(iter,fluid.mac_grid.grid){const VectorDi& cell=iter.Coord();
		//		if(//cell[0]==0||cell[0]==fluid.mac_grid.grid.cell_counts[0]-1||
		//			cell[1]==0||cell[1]==fluid.mac_grid.grid.cell_counts[1]-1)
		//			fluid.bc.Set_Psi_D(cell,(ushort)CellType::Solid);}
		//	fluid.bc.Set_Psi_D(VectorDi::Zero(),(ushort)CellType::Air);

		//}break;
		//case 5: {	////2d vortex sheet
		//	cell_counts = VectorDi::Ones() * s;
		//	fluid.Initialize(cell_counts, (real)length / cell_counts[0]);

		//	real v_scaler = (real)1.05;
		//	real radius = length/(real)4;
		//	VectorD center = fluid.mac_grid.grid.Center();

		//	iterate_face_in_one_dim(0, iter, fluid.mac_grid) {
		//		const VectorDi& face = iter.Coord();
		//		VectorD pos = fluid.mac_grid.Face_Center(0, face);
		//		if ((pos - center).norm() <= radius)fluid.velocity(0, face) = v_scaler * (center[1] - pos[1]);
		//	}

		//	iterate_face_in_one_dim(1, iter, fluid.mac_grid) {
		//		const VectorDi& face = iter.Coord();
		//		VectorD pos = fluid.mac_grid.Face_Center(1, face);
		//		if ((pos - center).norm() <= radius)fluid.velocity(1, face) = v_scaler * (pos[0] - center[0]);
		//	}
		//}break;
		case 6: {	////Taylor vortex
			json j_projection; j_projection = { {"update_A",false} };
			json j_fluid = { { "projection", j_projection } };
			fluid = std::make_shared<FluidEuler<d>>(j_fluid);
			real PI = atan(1) * (real)4;
			cell_counts[1] = cell_counts[0]; length = (real)2 * PI;
			fluid->Initialize(cell_counts, (real)length / cell_counts[0]);

			VectorD vortex_p1, vortex_p2;
			// two vortices are 0.81 apart
			vortex_p1[0] = fluid->mac_grid.grid.Center()[0] + (real)0.405;
			vortex_p2[0] = fluid->mac_grid.grid.Center()[0] - (real)0.405;
			vortex_p1[1] = fluid->mac_grid.grid.Center()[1]; vortex_p2[1] = vortex_p1[1];
			Set_Velocity_Taylor_Vortex(vortex_p1, vortex_p2);
		}break;
		//case 7: { //oblique plume
		//	fluid.use_maccormack = false;
		//	fluid.use_zero_extrapolation = false;
		//	fluid.Initialize(cell_counts, (real)length / cell_counts[0]);
		//	fluid.velocity.Fill((real)0.0, 0);

		//	VectorD source_center = fluid.mac_grid.grid.Position(AuxFunc::V<d>((real).1, (real).25, (real).17));
		//	real source_size = (real).025;
		//	Sphere<d>source(source_center, source_size);

		//	VectorD source_velocity = AuxFunc::V<d>((real)1, (real).32, (real).5);
		//	iterate_cell(iter, fluid.mac_grid.grid) {
		//		const VectorDi& cell = iter.Coord();
		//		const VectorD& pos = fluid.mac_grid.grid.Center(cell);
		//		if (source.Phi(pos) < 0) {
		//			fluid.bc.Set_Psi_D(cell, (ushort)CellType::Solid);
		//			for (int axis = 0; axis < d; axis++) {
		//				for (int side = 0; side < 2; side++) {
		//					VectorDi face = fluid.mac_grid.Cell_Incident_Face(axis, cell, side);
		//					fluid.bc.Set_Psi_N(axis, face, source_velocity[axis]);
		//				}
		//			}
		//		}
		//	}
		//}break;
		default: {Assert(false, "FluidEulerDriver not a valid case"); }break;
		}
		fluid->Enforce_Boundary_Conditions();
	}

	///////////////////////////////////////////////////////////////////////////
	////velocity setup
	void Set_Velocity_Taylor_Vortex(VectorD& vortex_p1, VectorD& vortex_p2) { //only works for 2d
		/*function [wz] = taylor_vortex(obj)
		pr1 = (obj.px-obj.sizex/2.-0.4).^2+(obj.py-obj.sizey/2.).^2;
		pr2 = (obj.px-obj.sizex/2.+0.4).^2+(obj.py-obj.sizey/2.).^2;
		wz = 1/0.3 * (2-pr1/0.09).*exp(0.5*(1-pr1/0.09));
		wz = wz+1/0.3 * (2-pr2/0.09).*exp(0.5*(1-pr2/0.09));
		end*/
		Field<real, d> wz;
		wz.Resize(fluid->mac_grid.grid.cell_counts);

		iterate_cell(iter, fluid->mac_grid.grid) {
			const VectorDi& cell = iter.Coord();
			const VectorD& pos = fluid->mac_grid.grid.Center(cell);
			real pr1 = pow((pos[0] - vortex_p1[0]), 2) + pow((pos[1] - vortex_p1[1]), 2);
			real pr2 = pow((pos[0] - vortex_p2[0]), 2) + pow((pos[1] - vortex_p2[1]), 2);
			wz(cell) = (real)1 / (real)0.3 * ((real)2 - pr1 / (real)0.09) * exp((real)0.5 * ((real)1 - pr1 / (real)0.09));
			wz(cell) += (real)1 / (real)0.3 * ((real)2 - pr2 / (real)0.09) * exp((real)0.5 * ((real)1 - pr2 / (real)0.09));
		}

		Poisson<d> poisson;
		poisson.Init(fluid->mac_grid,1000,1e-5);
		poisson.Update_Unknown([=](const VectorDi& cell)->bool {return true;});
		poisson.Update_Vol([&](const int axis,const VectorDi& face)->Scalar {
			if (fluid->mac_grid.Is_Axial_Boundary_Face(axis,face)) { return 0.; }
			else { return 1.; } 
			}); //can put dx here
		poisson.Send_To_Device();
		poisson.Update_b(wz);
		poisson.Solve();
		Field<real, d> sol = poisson.x;

		iterate_face(axis, iter, fluid->mac_grid) {
			const VectorDi& face = iter.Coord();
			if (fluid->mac_grid.Is_Boundary_Face(axis, face)) { continue; }
			const VectorDi cell_left = MacGrid<d>::Face_Left_Cell(axis, face);
			const VectorDi cell_right = MacGrid<d>::Face_Right_Cell(axis, face);
			if (axis == 0) { fluid->velocity(1, face) = (sol(cell_left)-sol(cell_right)) * fluid->mac_grid.grid.dx; }
			else if (axis == 1) { fluid->velocity(0, face) = (sol(cell_right)-sol(cell_left)) * fluid->mac_grid.grid.dx; }
		}
	}
};

#endif
