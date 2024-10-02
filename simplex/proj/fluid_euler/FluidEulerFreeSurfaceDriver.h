//////////////////////////////////////////////////////////////////////////
// Euler fluid with free surface driver
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __FluidEulerFreeSurfaceDriver_h__
#define __FluidEulerFreeSurfaceDriver_h__
#include "Common.h"
#include "File.h"
#include "GeometryPrimitives.h"
#include "FluidEulerFreeSurface.h"
#include "Particles.h"
#include "MarchingCubes.h"
#include "ImplicitShape.h"
#include "Driver.h"
#include "EulerInitializer.h"

template<int d> class FluidEulerFreeSurfaceDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;
public:
	std::shared_ptr<FluidEulerFreeSurface<d>> fluid;

	FluidEulerFreeSurfaceDriver() {}

	real CFL() const
	{
		real epsilon=(real)1e-5;
		return cfl * fluid->mac_grid.grid.dx / (fluid->Max_Abs(fluid->velocity) + epsilon);
	}

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		fluid->Advance(time, dt);
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
		//fluid->velocity.Write_To_File_3d(frame_dir + "/velocity");
		if (frame == 0) {
			////Write BC
			RenderFunc::Write_Dirichlet_Boundary_Conditions(frame_dir + "/psi_D", mac_grid, fluid->bc);
			RenderFunc::Write_Neumann_Boundary_Conditions(frame_dir + "/psi_N", mac_grid, fluid->bc);
		}

		////Write fluid type
		//{std::string file_name=frame_dir+"/fluid_type";
		//Field<real,d> fluid_type;fluid_type.Resize(mac_grid.grid.cell_counts);
		//int cell_num= mac_grid.grid.Number_Of_Cells();
		//#pragma omp parallel for
		//for(int i=0;i<cell_num;i++){const VectorDi& cell= mac_grid.grid.Cell_Coord(i);
		//	if(fluid->type(cell)==(ushort)CellType::Fluid)fluid_type(cell)=(real)0;
		//	else if(fluid->type(cell)==(ushort)CellType::NB)fluid_type(cell)=(real).25;
		//	else if(fluid->type(cell)==(ushort)CellType::BD)fluid_type(cell)=(real).5;
		//	else if(fluid->type(cell)==(ushort)CellType::Air)fluid_type(cell)=(real).75;
		//	else fluid_type(cell)=(real)1;}
		//fluid_type.Write_To_File_3d(file_name);}

		////Write interface
		//fluid->levelset.phi.Write_To_File_3d(frame_dir + "/phi");

		//if (frame != 0)fluid->projection->cell_b.Write_To_File_3d(frame_dir + "/solid_phi");
		
		RenderFunc::Write_LevelSet_Mesh(frame_dir + (d == 2 ? "/segment_mesh" : "/triangle_mesh"), fluid->levelset);

		if (fluid->spray_system) {
			RenderFunc::Write_Points_Float<d, real>(frame_dir + "/trial_points", fluid->spray_system->trial_points.XRef());
			RenderFunc::Write_Points_Float<d, real>(frame_dir + "/spray_points", fluid->spray_system->spray_points.XRef());
			RenderFunc::Write_Points_Float<d, real>(frame_dir + "/foam_points", fluid->spray_system->foam_points.XRef());
		}

		std::cout<<"Write to files, frame "<<frame<<std::endl;
		AuxFunc::Seperation_Line();
	}
	
	virtual void Initialize()
	{
		int s=scale;real length=(real)1;VectorDi cell_counts=VectorDi::Ones();cell_counts[0]=s;
		cell_counts[1]=cell_counts[0]/2;if(d>2)cell_counts[2]=cell_counts[0]/2;

		switch(test){
		case 1:{	////droplet falling into water tank
			VectorD g_coeff = VectorD::Unit(1) * (-9.8);
			json j_projection = {
				{"surface_tension_mode",SurfaceTensionMode::IMPLICIT},
				{"sigma",0.01},
				{"use_jump_gravity",true},
				{"gravity_coefficient",g_coeff},
			};
			json j_fluid = {
				{"projection", j_projection},
				{"verbose",false},
				{"use_body_force",false},
				{"narrow_band_cell_num",3},
				{"dirac_band_cell_num",3}
			};
			fluid = std::make_shared<FluidEulerFreeSurface<d>>(j_fluid);
			EulerInitializer<d> perimeter;
			perimeter.Set_Domain(length, scale, VectorDi::Ones());
			perimeter.Set_Boundary_Width(0, 0, 0, 0, 0, 0);
			perimeter.Set_Boundary_Value(0, 0, 0, 0, 0, 0);
			perimeter.Generate_Parameters();
			fluid->Initialize(perimeter.cell_counts, perimeter.dx, perimeter.domain_min);
			perimeter.Fill_Boundary_Condition(fluid->bc);

			//Initialize levelset
			VectorD center = AuxFunc::V<d>(0, 0.3, 0) * length;
			auto sphere_ptr = std::make_shared<Sphere<d>>(center, 0.1 * length);
			auto plane_ptr = std::make_shared<Plane<d>>(VectorD::Unit(1), AuxFunc::V<d>(0, 0, 0) * length);
			ImplicitShape<d> water;
			water += sphere_ptr;
			water += plane_ptr;
			fluid->levelset.Set_By_Shape(water);
		}break;
		case 2:{	////oscillating droplet with surface tension
			////use_jump_condition true:explicit surface tension/magnetic; false:implicit surface tension/magnetic
			// Issue: .\fluid_euler.exe -driver 2 -test 2 -o output -s 192 -lf 20 -d 3
			json j_projection = { {"surface_tension_mode",SurfaceTensionMode::IMPLICIT}, {"sigma",4e4} };
			json j_fluid = {
				{ "projection", j_projection } ,
				{"verbose",false},
				{"use_body_force",false},
				{"narrow_band_cell_num",3},
				{"dirac_band_cell_num",3}
			};
			length = (real)100;
			fluid = std::make_shared<FluidEulerFreeSurface<d>>(j_fluid);
			cell_counts=VectorDi::Ones()*scale;
			fluid->Initialize(cell_counts,(real)length/cell_counts[0]);

			const auto mac_grid = fluid->mac_grid;
			VectorD e_r = VectorD::Ones() * (real)8;
			e_r[1] *= 2;
			Ellipsoid<d> ellipsoid(mac_grid.grid.Center(),e_r);
			fluid->levelset.Set_By_Geom(ellipsoid);
			fluid->levelset.Fast_Marching();
		}break;
		case 3:{	////droplet falling onto ground
			json j_projection = { {"surface_tension_mode",SurfaceTensionMode::EXPLICIT} ,{"sigma",1.0},{"use_vol_control",true} };
			json j_fluid = { { "projection", j_projection } ,{"verbose",false}, {"use_body_force",true} };
			fluid = std::make_shared<FluidEulerFreeSurface<d>>(j_fluid);
			cell_counts=VectorDi::Ones()*scale;
			fluid->Initialize(cell_counts,(real)length/cell_counts[0]);
			//fluid.use_body_force=true;
			//fluid.projection.use_surface_tension = true;
			//fluid.projection.use_jump_condition=false;	////true:explicit surface tension/magnetic; false:implicit surface tension/magnetic
			//fluid.projection.sigma=(real)1;
			frame_rate=200;

			const auto mac_grid = fluid->mac_grid;
			VectorD e_r=VectorD::Ones()* mac_grid.grid.Length()[1]*(real).1;
			VectorD pos= mac_grid.grid.Center();pos[1]*=(real).5;
			Ellipsoid<d> ellipsoid(pos,e_r);
			iterate_cell(iter, mac_grid.grid){const VectorDi& cell=iter.Coord();
				fluid->levelset.phi(cell)=ellipsoid/*sphere*/.Phi(mac_grid.grid.Center(cell));}
			fluid->levelset.Fast_Marching();

			////initialize droplet velocity
			VectorD init_droplet_vel=VectorD::Unit(1)*(real)-2;
			fluid->velocity.Fill((real)0);
			iterate_face(axis, iter, mac_grid) {
				const VectorDi& face = iter.Coord();
				VectorD pos = mac_grid.Face_Center(axis, face);
				real phi = ellipsoid.Phi(pos);
				if (phi < fluid->narrow_band_width) {
					fluid->velocity(axis, face) = init_droplet_vel[axis];
				}
			}

			//fluid.projection.use_vol_control=true;
			auto free_proj = std::dynamic_pointer_cast<ProjectionIrregularCPX<real, d>>(fluid->projection);
			free_proj->target_vol = free_proj->current_vol = fluid->levelset.Total_Volume();

			////wall boundary
			for (int axis = 0; axis < d; axis++)iterate_face_in_one_dim(axis, iter, mac_grid) {
				const VectorDi& face = iter.Coord();
				if (face[axis] == 0 || face[axis] == mac_grid.face_grids[axis].node_counts[axis] - 1)fluid->bc.Set_Psi_N(axis, face, (real)0);
			}
		}break;
		case 4: 
  		{ // karman vortex in sink with source
			Assert(d == 3, "Case 4 only available in 3D");
			length = 1;
			cell_counts = VectorDi::Ones() * scale;
			cell_counts[1] /= 4;
			cell_counts[2] /= 2;
			real length0 = length * (real)cell_counts[0]/(real)cell_counts[0];
			real length1 = length * (real)cell_counts[1]/(real)cell_counts[0];
			real length2 = length * (real)cell_counts[2]/(real)cell_counts[0];

			real dx = (real)length / cell_counts[0];
			json j_projection = {{"surface_tension_mode",SurfaceTensionMode::NONE}};
			json j_fluid = {{"projection",j_projection},{"use_body_force",true},{"sigma",(real)0.}};
			fluid = std::make_shared<FluidEulerFreeSurface<d>>(j_fluid);
			fluid->Initialize(cell_counts, dx);
			const auto mac_grid = fluid->mac_grid;
			
		

			

			////Source
			real source_speed = (real)1; int axis = 0;
			iterate_face_in_one_dim(axis, iter, mac_grid) {
				const VectorDi& face = iter.Coord();
				if (face[0] < 5)fluid->bc.Set_Psi_N(0, face, source_speed);
			}

			////Wall
			{axis = 1; iterate_face_in_one_dim(axis, iter, mac_grid) {
				const VectorDi& face = iter.Coord();
				if (face[axis] == 0 )fluid->bc.Set_Psi_N(axis, face, 0);
			}}
			if (d == 3) {
				axis = 2; iterate_face_in_one_dim(axis, iter, mac_grid) {
					const VectorDi& face = iter.Coord();
					if (face[axis] == 0 || face[axis] == mac_grid.face_grids[axis].node_counts[axis] - 1)fluid->bc.Set_Psi_N(axis, face, 0);
				}
			}
					 


			frame_rate = 500;
			////Water bulk
			
			Plane<d> plane(VectorD::Unit(1), mac_grid.grid.Center() );
			fluid->levelset.Set_By_Geom(plane);
			fluid->levelset.Fast_Marching();


			//Solid
			VectorD solid_center = mac_grid.grid.Center();
			solid_center[0] = length * 0.2;
			solid_center[1] = length1 * 0.75;
			VectorD tube_normal =  VectorD::Unit(1);
			real tube_radius = 0.04 * length;
			real tube_height = 0.6 * length1;
			Capsule<d> obstacle(solid_center, tube_radius, tube_normal, tube_height);
			iterate_cell(iter, mac_grid.grid) {
				const VectorDi& cell = iter.Coord(); const VectorD& pos = mac_grid.grid.Center(cell);
				if (obstacle.Inside(pos)) {
					fluid->bc.Set_Psi_D(cell, (ushort)CellType::Solid);
					for (int axis = 0; axis < d; axis++)for (int side = 0; side < 2; side++) { VectorDi face = mac_grid.Cell_Incident_Face(axis, cell, side); fluid->bc.Set_Psi_N(axis, face, (real)0); }
				}
			}

			fluid->Update_Cell_Types();
		}break;
		case 5: { // Water ripples in a free surface flow
			length = 1;
			cell_counts = VectorDi::Ones() * scale;
			real length0 = length * cell_counts[0] / cell_counts[0];
			real length1 = length * cell_counts[1] / cell_counts[0];
			real length2 = length * cell_counts[2] / cell_counts[0];
			real dx = length / cell_counts[0];
			real sigma = 0e-1;
			VectorD gcoeff = VectorD::Unit(1) * -9.8;
			json j_projection = {
				{"surface_tension_mode",SurfaceTensionMode::EXPLICIT},
				{"sigma",sigma},
				{"use_jump_gravity",true},
				{"gravity_coefficient",gcoeff}
			};
			json j_fluid = {
				{"projection",j_projection},
				{"h_bar",0.2},
				{"use_solid",false},
				{"use_body_force",false},
				{"use_implicit_surface_tension",false}
			};
			fluid = std::make_shared<FluidEulerFreeSurface<d>>(j_fluid);
			fluid->Initialize(cell_counts, (real)length / cell_counts[0]);

			EulerInitializer<d> perimeter;
			perimeter.Set_Boundary_Width(0, 0, 0, -1, 0, 0);
			perimeter.Set_Boundary_Value(0, 0, 0, 0, 0, 0);
			perimeter.Set_Parameters(fluid->mac_grid);
			perimeter.Fill_Boundary_Condition(fluid->bc);

			frame_rate = 200;

			real ly = 0.8 * length;
			int ny = floor(ly/dx);
			int nyh = floor(ny/1.2);
			int nyl = floor(ny/2);
			ly = ny * dx;

			for (int axis = 0; axis < d; axis++) iterate_face_in_one_dim(axis, iter, fluid->mac_grid) {
				const VectorDi& face = iter.Coord();
				VectorD pos = fluid->mac_grid.Face_Center(axis, face);
				VectorD out_center = fluid->mac_grid.grid.Center();
				out_center[1] = pos[1];
				real out_r = length * 0.04;
				bool is_out = (pos - out_center).squaredNorm() > out_r * out_r;
				if (is_out && face[1] == ny && axis == 1) fluid->bc.Set_Psi_N(axis, face, (real)0);
				if (axis != 1 && face[1] < nyh && face[1] > nyl) fluid->bc.Del_Psi_N(axis,face);
			}

			////Water bulk
			VectorD plane_center1 = fluid->mac_grid.grid.Center();
			plane_center1[1] = ly;
			Plane<d> plane1(VectorD::Unit(1), plane_center1);

			VectorD plane_center2 = fluid->mac_grid.grid.Center();
			plane_center2[1] = length * 0.98;
			Plane<d> plane2(VectorD::Unit(1), plane_center2);

			VectorD plane_center3 = fluid->mac_grid.grid.Center();
			plane_center3[1] = length * 0.3;
			Plane<d> plane3(VectorD::Unit(1), plane_center3);

			iterate_cell(iter, fluid->mac_grid.grid) {
				const VectorDi& cell = iter.Coord();
				const VectorD& pos = fluid->mac_grid.grid.Center(cell);
				fluid->levelset.phi(cell) = std::min(plane1.Phi(pos) * plane2.Phi(pos),plane3.Phi(pos));
			}

			fluid->levelset.Fast_Marching();

			fluid->Update_Cell_Types();
		}break;
		case 6: {//droplet falling with spray system
			EulerInitializer<d> perimeter;
			perimeter.Set_Domain(length, scale, VectorDi::Ones());
			perimeter.Set_Boundary_Width(0, 0, 0, 0, 0, 0);
			perimeter.Set_Boundary_Value(0, 0, 0, 0, 0, 0);
			perimeter.Generate_Parameters();

			VectorD g_coeff = VectorD::Unit(1) * (-9.8);
			json j_projection = {
				{"surface_tension_mode",SurfaceTensionMode::NONE},
				{"sigma",0.01},
				{"use_jump_gravity",true},
				{"gravity_coefficient",g_coeff},
			};
			json j_spray = {
				{"phi_gen",-1.5 * perimeter.dx},//note it's a negative number
				{"cell_sample_num",1},
				{"gravity_acc",g_coeff},
				{"turn_form_prob",0.5},
				{"trial_die_time",1.0},
				{"foam_die_time",1.0},
			};
			json j_fluid = {
				{"projection", j_projection},
				{"spray_system",j_spray},
				{"verbose",false},
				{"use_body_force",false},
				{"narrow_band_cell_num",3},
				{"dirac_band_cell_num",3}
			};
			fluid = std::make_shared<FluidEulerFreeSurface<d>>(j_fluid);
			fluid->Initialize(perimeter.cell_counts, perimeter.dx, perimeter.domain_min);
			perimeter.Fill_Boundary_Condition(fluid->bc);

			//Initialize levelset
			VectorD center = AuxFunc::V<d>(0, 0.3, 0) * length;
			auto sphere_ptr = std::make_shared<Sphere<d>>(center, 0.1 * length);
			auto plane_ptr = std::make_shared<Plane<d>>(VectorD::Unit(1), AuxFunc::V<d>(0, 0, 0) * length);
			ImplicitShape<d> water;
			water += sphere_ptr;
			water += plane_ptr;
			fluid->levelset.Set_By_Shape(water);
		}break;
		case 7: {//dam break
			EulerInitializer<d> perimeter;
			perimeter.Set_Domain(length, scale, VectorDi::Ones());
			perimeter.Set_Boundary_Width(0, 0, 0, 0, 0, 0);
			perimeter.Set_Boundary_Value(0, 0, 0, 0, 0, 0);
			perimeter.Generate_Parameters();

			VectorD g_coeff = VectorD::Unit(1) * (-9.8);
			json j_projection = {
				{"surface_tension_mode",SurfaceTensionMode::NONE},
				{"sigma",0.01},
				{"use_jump_gravity",true},
				{"gravity_coefficient",g_coeff},
			};
			json j_spray = {
				{"phi_gen",-1.5 * perimeter.dx},//note it's a negative number
				{"cell_sample_num",1},
				{"gravity_acc",g_coeff},
				{"turn_form_prob",0.5},
				{"trial_die_time",1.0},
				{"foam_die_time",1.0},
			};
			json j_fluid = {
				{"projection", j_projection},
				{"spray_system",j_spray},
				{"verbose",false},
				{"use_body_force",false},
				{"narrow_band_cell_num",3},
				{"dirac_band_cell_num",3}
			};
			fluid = std::make_shared<FluidEulerFreeSurface<d>>(j_fluid);
			fluid->Initialize(perimeter.cell_counts, perimeter.dx, perimeter.domain_min);
			perimeter.Fill_Boundary_Condition(fluid->bc);

			//Initialize levelset
			VectorD min_crn = AuxFunc::V<d>(-0.5, -0.5, -0.5) * length;
			VectorD max_crn = AuxFunc::V<d>(-0.25, 0, 0) * length;
			Box<d> box(min_crn, max_crn);
			fluid->levelset.Set_By_Geom(box);
		}break;
		default: {Assert(false, "FluidEulerDriver not a valid case"); }break;
		}
	}
};

#endif