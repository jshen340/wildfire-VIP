//////////////////////////////////////////////////////////////////////////
// Flip Fluid driver
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __FluidPicFlipDriver_h__
#define __FluidPicFlipDriver_h__
#include "Common.h"
#include "File.h"
#include "Timer.h"
#include "GeometryPrimitives.h"
#include "FluidPicFlipOld.h"
#include "Particles.h"
#include "Driver.h"
#include "RenderFunc.h"

template<int d> class FluidPicFlipDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;
public:
	FluidPicFlipOld<d> fluid;

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		fluid.Advance(dt);
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Timer timer;
		timer.Reset();
		Base::Write_Output_Files(frame);
		if(frame==0){
			std::string file_name=frame_dir+"/grid";
			fluid.mac_grid.grid.Write_To_File_3d(file_name);
			std::cout<<"Write to file "<<file_name<<std::endl;}

		if(!File::Directory_Exists(frame_dir.c_str()))File::Create_Directory(frame_dir);
		
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
		particles.Write_To_File_3d(file_name);}
		
		////Write fluid type
		{std::string file_name=frame_dir+"/fluid_type";
		Field<real,d> fluid_type;fluid_type.Resize(fluid.mac_grid.grid.cell_counts);
		iterate_cell(iter,fluid.mac_grid.grid){const VectorDi& cell=iter.Coord();
			if(fluid.type(cell)==(ushort)CellType::Fluid)fluid_type(cell)=(real)0;
			else fluid_type(cell)=(real)1;}
		fluid_type.Write_To_File_3d(file_name);}

		////Write Particles
		RenderFunc::Write_Points_Float<d, real>(frame_dir + "/tracker_points", fluid.points.XRef());
		timer.Elapse_And_Output("IO");
	}
	
	virtual void Initialize()
	{
		int s=scale;real length=(real)1;VectorDi cell_counts=VectorDi::Ones();cell_counts[0]=s;
		cell_counts[1]=cell_counts[0]/2;if(d>2)cell_counts[2]=cell_counts[0]/2;

		switch(test){
		case 1:{	////droplet falling into water
			cell_counts=VectorDi::Ones()*scale;
			fluid.Initialize(cell_counts,(real)length/cell_counts[0]);

			//fluid.use_interface=true;
			//fluid.Initialize_Interface();
			fluid.use_body_force=true;
			
			////droplet
			VectorD center=fluid.mac_grid.grid.Center();
			center[1]=fluid.mac_grid.grid.domain_min[1]+fluid.mac_grid.grid.Length()[1]*(real).8;
			real r=fluid.mac_grid.grid.Length()[1]*(real).1;
			Sphere<d> sphere(center,r);
			////water bulk
			Plane<d> plane(VectorD::Unit(1),fluid.mac_grid.grid.Center());
			////initialize phi
			iterate_cell(iter,fluid.mac_grid.grid){const VectorDi& cell=iter.Coord();
				VectorD pos=fluid.mac_grid.grid.Center(cell);
				real phi=std::min(plane.Phi(pos),sphere.Phi(pos));
				if(phi<(real)0)fluid.Seed_Particles(cell,(int)pow(2,d));}

			////wall boundary
			for(int axis=0;axis<d;axis++)iterate_face_in_one_dim(axis,iter,fluid.mac_grid){const VectorDi& face=iter.Coord();
				if(face[axis]==0||face[axis]==fluid.mac_grid.face_grids[axis].node_counts[axis]-1)fluid.bc.Set_Psi_N(axis,face,(real)0);}

			//fluid.projection.solver_mode = SolverType::KRYLOV_CPU;
		}break;}
		fluid.Enforce_Boundary_Conditions();
	}
};

#endif
