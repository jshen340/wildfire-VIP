//////////////////////////////////////////////////////////////////////////
// Euler fluid with two phase flow driver
// Copyright (c) (2018-), Bo Zhu, Shiying Xiong, Zhecheng Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __FluidEulerTwoPhaseDriver_h__
#define __FluidEulerTwoPhaseDriver_h__
#include "Common.h"
#include "File.h"
#include "GeometryPrimitives.h"
#include "FluidEulerTwoPhase.h"
#include "Particles.h"
#include "MarchingCubes.h"
#include "Driver.h"
#include "EulerInitializer.h"

template<int d> class FluidEulerTwoPhaseDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;
public:
	std::shared_ptr<FluidEulerTwoPhase<d>> fluid;

	FluidEulerTwoPhaseDriver() {}

	real CFL() const
	{
		real epsilon=(real)1e-5;
		return cfl * fluid->mac_grid.grid.dx / (fluid->Max_Abs(fluid->velocity) + epsilon);
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
		if (frame == 0) {
			////Write BC
			RenderFunc::Write_Dirichlet_Boundary_Conditions(frame_dir + "/psi_D", mac_grid, fluid->bc);
			RenderFunc::Write_Neumann_Boundary_Conditions(frame_dir + "/psi_N", mac_grid, fluid->bc);
		}

		////Write fluid type
		{std::string file_name=frame_dir+"/fluid_type";
		Field<real,d> fluid_type;fluid_type.Resize(mac_grid.grid.cell_counts);
		int cell_num= mac_grid.grid.Number_Of_Cells();
		#pragma omp parallel for
		for(int i=0;i<cell_num;i++){const VectorDi& cell= mac_grid.grid.Cell_Coord(i);
			if(fluid->type(cell)==(ushort)CellType::Liquid)fluid_type(cell)=(real)0;
			else if(fluid->type(cell)==(ushort)CellType::NB)fluid_type(cell)=(real).25;
			else if(fluid->type(cell)==(ushort)CellType::BD)fluid_type(cell)=(real).5;
			else if(fluid->type(cell)==(ushort)CellType::Air)fluid_type(cell)=(real).75;
			else fluid_type(cell)=(real)1;}
		fluid_type.Write_To_File_3d(file_name);}

		////Write interface
		fluid->levelset.phi.Write_To_File_3d(frame_dir + "/phi");
		
		RenderFunc::Write_LevelSet_Mesh(frame_dir + (d == 2 ? "/segment_mesh" : "/triangle_mesh"), fluid->levelset);

		std::cout<<"Write to files, frame "<<frame<<std::endl;
		AuxFunc::Seperation_Line();
	}
	
	virtual void Initialize()
	{
		int s=scale;real length=(real)1;VectorDi cell_counts=VectorDi::Ones();cell_counts[0]=s;
		cell_counts[1]=cell_counts[0]/2;if(d>2)cell_counts[2]=cell_counts[0]/2;

		switch(test){
		case 1:{	////droplet falling into water tank
			////use_jump_condition true:explicit surface tension/magnetic; false:implicit surface tension/magnetic
			VectorD g_coeff = VectorD::Unit(1) * (-9.8);
			json j_projection = {
				{"surface_tension_mode",SurfaceTensionMode::EXPLICIT},
				{"sigma",0.01},
				{"use_jump_gravity",true},
				{"gravity_coefficient",g_coeff}
			};
			json j_fluid = { 
				{"projection", j_projection}, 
				{"verbose",false}, 
				{"use_body_force",false} 
			};
			fluid = std::make_shared<FluidEulerTwoPhase<d>>(j_fluid);
			cell_counts=VectorDi::Ones()*scale;
			fluid->Initialize(cell_counts,(real)length/cell_counts[0]);

			////wall boundary
			EulerInitializer<d> perimeter;
			perimeter.Set_Boundary_Width(0, 0, 0, 0, 0, 0);
			perimeter.Set_Boundary_Value(0, 0, 0, 0, 0, 0);
			perimeter.Set_Parameters(fluid->mac_grid);
			perimeter.Fill_Boundary_Condition(fluid->bc);

			const auto mac_grid = fluid->mac_grid;
			////droplet
			VectorD center= mac_grid.grid.Center();
			center[1]= mac_grid.grid.domain_min[1]+ mac_grid.grid.Length()[1]*(real).8;
			real r= mac_grid.grid.Length()[1]*(real).1;
			Sphere<d> sphere(center,r);

			////water bulk
			Plane<d> plane(VectorD::Unit(1), mac_grid.grid.Center());
			
			////combined phi 
			iterate_cell(iter, mac_grid.grid) {
				const VectorDi& cell = iter.Coord();
				VectorD pos = mac_grid.grid.Center(cell);
				fluid->levelset.phi(cell) = std::min(plane.Phi(pos), sphere.Phi(pos));
			}
		}break;
		case 2:{	////oscillating droplet with surface tension
			json j_projection = { 
				{"surface_tension_mode",SurfaceTensionMode::EXPLICIT}, 
				{"sigma",1e-2}, 
			};
			json j_fluid = { 
				{"projection", j_projection}, 
				{"verbose",false}, 
				{"use_body_force",false} 
			};
			fluid = std::make_shared<FluidEulerTwoPhase<d>>(j_fluid);
			cell_counts=VectorDi::Ones()*scale;
			fluid->Initialize(cell_counts,(real)length/cell_counts[0]);

			const auto mac_grid = fluid->mac_grid;
			VectorD e_r=VectorD::Ones()* mac_grid.grid.Length()[1]*(real).2;//e_r[1]*=(real).5;
			Ellipsoid<d> ellipsoid(mac_grid.grid.Center(),e_r);
			fluid->levelset.Set_By_Geom(ellipsoid);
			fluid->levelset.Fast_Marching();
		}break;
		case 3:{	////droplet falling onto ground
			VectorD g_coeff = VectorD::Unit(1) * (-9.8);
			json j_projection = { 
				{"surface_tension_mode",SurfaceTensionMode::IMPLICIT},
				{"sigma",1.0}, 
				{"use_vol_control",false},
				{"use_jump_gravity",true},
				{"gravity_coefficient", g_coeff}
			};
			json j_fluid = { {"projection", j_projection}, {"verbose",false}, {"use_body_force",false} };
			fluid = std::make_shared<FluidEulerTwoPhase<d>>(j_fluid);
			cell_counts=VectorDi::Ones()*scale;
			fluid->Initialize(cell_counts,(real)length/cell_counts[0]);
			frame_rate=200;

			////wall boundary
			EulerInitializer<d> perimeter;
			perimeter.Set_Boundary_Width(0, 0, 0, 0, 0, 0);
			perimeter.Set_Boundary_Value(0, 0, 0, 0, 0, 0);
			perimeter.Set_Parameters(fluid->mac_grid);
			perimeter.Fill_Boundary_Condition(fluid->bc);

			const auto mac_grid = fluid->mac_grid;
			VectorD e_r=VectorD::Ones()* mac_grid.grid.Length()[1]*(real).1;
			VectorD pos= mac_grid.grid.Center();pos[1]*=(real).5;
			Ellipsoid<d> ellipsoid(pos,e_r);
			fluid->levelset.Set_By_Geom(ellipsoid);
			fluid->levelset.Fast_Marching();

			auto proj = std::dynamic_pointer_cast<ProjectionTwoPhaseCPX<real, d>>(fluid->projection);
			proj->target_vol = proj->current_vol = fluid->levelset.Total_Volume();
		}break;
		case 4:{	////ascneding water bubble in higher density fluid
			VectorD g_coeff = VectorD::Unit(1) * (-9.8);
			json j_projection = { 
				{"surface_tension_mode",SurfaceTensionMode::IMPLICIT}, 
				{"rhoA",2.}, 
				{"rhoL",1.}, 
				{"sigma",0e-4},
				{"use_jump_gravity",true},
				{"gravity_coeff",g_coeff}
			};
			json j_fluid = { 
				{"projection", j_projection}, 
				{"verbose",false}, 
				{"use_body_force",false}, 
				{"rhoA",2.}, 
				{"rhoL",1.} 
			};
			fluid = std::make_shared<FluidEulerTwoPhase<d>>(j_fluid);
			cell_counts=VectorDi::Ones()*scale;
			fluid->Initialize(cell_counts,(real)length/cell_counts[0]);

			////wall boundary
			EulerInitializer<d> perimeter;
			perimeter.Set_Boundary_Width(0, 0, 0, 0, 0, 0);
			perimeter.Set_Boundary_Value(0, 0, 0, 0, 0, 0);
			perimeter.Set_Parameters(fluid->mac_grid);
			perimeter.Fill_Boundary_Condition(fluid->bc);

			const auto mac_grid = fluid->mac_grid;

			////init levelset
			real r = length * (real).2;
			Sphere<d> sphere(mac_grid.grid.Center(), r);
			fluid->levelset.Set_By_Geom(sphere);
			fluid->levelset.Fast_Marching();
		}break;
		default: {Assert(false, "FluidEulerDriver not a valid case"); }break;
		}
	}
};

#endif