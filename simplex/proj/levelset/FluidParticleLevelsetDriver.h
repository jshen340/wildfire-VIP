#pragma once
#include "Common.h"
#include "File.h"
#include "GeometryPrimitives.h"
#include "FluidParticleLevelset.h"
#include "Particles.h"
#include "MarchingCubes.h"
#include "Driver.h"

template<int d> class FluidParticleLevelsetDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;
public:
	FluidParticleLevelset<d> fluid;

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		fluid.Advance(dt);
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);
		if(frame==0){
			std::string file_name=frame_dir+"/grid";
			fluid.mac_grid.grid.Write_To_File_3d(file_name);
			std::cout<<"Write to file "<<file_name<<std::endl;}
		
		////Write velocity
		{std::string file_name=frame_dir+"/velocity";
		fluid.velocity.Write_To_File_3d(file_name);}

		////Write BC
		{std::string file_name=frame_dir+"/psi_D";
		Particles<d> particles;
		for(auto p:fluid.bc.psi_D_values){
			VectorDi cell=fluid.mac_grid.grid.Cell_Coord(p.first);
			VectorD pos=fluid.mac_grid.grid.Center(cell);
			int i=particles.Add_Element();particles.X(i)=pos;}
		particles.Write_To_File_3d(file_name);}
		{std::string file_name=frame_dir+"/psi_N";
		Particles<d> particles;
		for(auto p:fluid.bc.psi_N_values){int axis=p.first[0];
			VectorDi face=fluid.mac_grid.Face_Coord(axis,p.first[1]);
			VectorD pos=fluid.mac_grid.Face_Center(axis,face);
			int i=particles.Add_Element();particles.X(i)=pos;}
		File::Write_Binary_To_File(file_name,particles);
		particles.Write_To_File_3d(file_name);}

		////Write fluid type
		{std::string file_name=frame_dir+"/fluid_type";
		Field<real,d> fluid_type;fluid_type.Resize(fluid.mac_grid.grid.cell_counts);
		iterate_cell(iter,fluid.mac_grid.grid){const VectorDi& cell=iter.Coord();
			if(fluid.type(cell)==(ushort)CellType::Fluid)fluid_type(cell)=(real)0;
			else fluid_type(cell)=(real)1;}
		fluid_type.Write_To_File_3d(file_name);}

		////Write interface
		std::string file_name=frame_dir+"/phi";
		fluid.levelset.phi.Write_To_File_3d(file_name);

		MarchingCubes<d> marching_cubes(fluid.levelset);
		marching_cubes.Marching();
		if(d==2){
			std::string file_name=frame_dir+"/segment_mesh";
			(*marching_cubes.mesh).Write_To_File_3d(file_name);}
		else{
			std::string file_name=frame_dir+"/triangle_mesh";
			(*marching_cubes.mesh).Write_To_File_3d(file_name);}

        {std::string file_name=frame_dir+"/particles";
		    fluid.pls.particles.Write_To_File_3d(file_name);}

		{std::string file_name = frame_dir + "/point_impulse";
		for(int i=0;i<fluid.pls.particles.Size();i++){
			VectorD norm=fluid.pls.Normal(fluid.pls.particles.X(i));
			real r=fluid.pls.Particle_Radius(i);
			fluid.pls.particles.F(i)=-fluid.pls.Particle_Sign(i)*norm.normalized()*r;}
		Write_Segments_To_File_3d_Fast<d, real>(fluid.pls.particles.XRef(), fluid.pls.particles.FRef(), file_name);
		}

		{std::string file_name = frame_dir + "/point_velocity";
		Write_Segments_To_File_3d_Fast<d, real>(fluid.pls.particles.XRef(), fluid.pls.particles.VRef(), file_name);
		}
		
	}
	
	virtual void Initialize()
	{
		int s=scale;real length=(real)1;VectorDi cell_counts=VectorDi::Ones();cell_counts[0]=s;
		cell_counts[1]=cell_counts[0]/2;if(d>2)cell_counts[2]=cell_counts[0]/2;

		switch(test){
		case 1:{	////Droplet falling into water
			cell_counts=VectorDi::Ones()*scale;
			fluid.Initialize(cell_counts,(real)length/cell_counts[0]);
			fluid.use_body_force=true;
			
			////Droplet
			VectorD center=fluid.mac_grid.grid.Center();
			center[1]=fluid.mac_grid.grid.domain_min[1]+fluid.mac_grid.grid.Length()[1]*(real).8;
			real r=fluid.mac_grid.grid.Length()[1]*(real).1;
			Sphere<d> sphere(center,r);
			////Water bulk
			Plane<d> plane(VectorD::Unit(1),fluid.mac_grid.grid.Center());
			////Initialize phi
			iterate_cell(iter,fluid.mac_grid.grid){const VectorDi& cell=iter.Coord();
				VectorD pos=fluid.mac_grid.grid.Center(cell);
				fluid.levelset.phi(cell)=std::min(plane.Phi(pos),sphere.Phi(pos));}
            
            fluid.Setup_Particle_Level_Set();

			////Wall boundary
			for(int axis=0;axis<d;axis++)iterate_face_in_one_dim(axis,iter,fluid.mac_grid){const VectorDi& face=iter.Coord();
				if(face[axis]==0||face[axis]==fluid.mac_grid.face_grids[axis].node_counts[axis]-1)fluid.bc.Set_Psi_N(axis,face,(real)0);}

			fluid.projection.use_multigrid_solver=true;
		}break;
		case 2:{	////Surface tension
			cell_counts=VectorDi::Ones()*scale;
			fluid.Initialize(cell_counts,(real)length/cell_counts[0]);
			fluid.use_body_force=false;
			fluid.projection.use_jump_condition=true;
            fluid.projection.sigma=5e-3;

			VectorD e_r=VectorD::Ones()*fluid.mac_grid.grid.Length()[1]*(real).2;e_r[1]*=(real).5;
			Ellipsoid<d> ellipsoid(fluid.mac_grid.grid.Center(),e_r);
			//Sphere<d> sphere(fluid.mac_grid.grid.Center(),e_r[0]);
			iterate_cell(iter,fluid.mac_grid.grid){const VectorDi& cell=iter.Coord();
				fluid.levelset.phi(cell)=ellipsoid/*sphere*/.Phi(fluid.mac_grid.grid.Center(cell));}
			fluid.levelset.Fast_Marching();
			fluid.Setup_Particle_Level_Set();
		}break;}

		fluid.Enforce_Boundary_Conditions();
		fluid.Advance((real).01);
	}

	
};