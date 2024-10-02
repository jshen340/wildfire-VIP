//////////////////////////////////////////////////////////////////////////
// Euler fluid immersed boundary driver
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __FluidEulerImmersedBoundaryDriver_h__
#define __FluidEulerImmersedBoundaryDriver_h__
#include "Common.h"
#include "File.h"
#include "GeometryPrimitives.h"
#include "FluidEulerImmersedBoundaryOld.h"
#include "Particles.h"
#include "MarchingCubes.h"
#include "Driver.h"

template<int d> class FluidEulerImmersedBoundaryDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;
public:
	FluidEulerImmersedBoundaryOld<d> fluid;

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
	}
	
	virtual void Initialize()
	{
		int s=scale;real length=(real)1;VectorDi cell_counts=VectorDi::Ones();cell_counts[0]=s;
		cell_counts[1]=cell_counts[0]/2;if(d>2)cell_counts[2]=cell_counts[0]/2;
		
		switch(test){
		case 1:{
			fluid.Initialize(cell_counts,(real)length/cell_counts[0]);
			fluid.velocity.Fill((real).2,0);

			////Source
			real source_speed=(real)1;int axis=0;
			iterate_face_in_one_dim(axis,iter,fluid.mac_grid){const VectorDi& face=iter.Coord();
				if(face[0]==0||face[0]==fluid.mac_grid.face_grids[0].node_counts[0]-1)fluid.bc.Set_Psi_N(0,face,source_speed);}
			////Wall
			axis=1;iterate_face_in_one_dim(axis,iter,fluid.mac_grid){const VectorDi& face=iter.Coord();
				if(face[axis]==0||face[axis]==fluid.mac_grid.face_grids[axis].node_counts[axis]-1)fluid.bc.Set_Psi_N(axis,face,0);}
			////Solid

			Sphere<d> sphere(fluid.mac_grid.grid.Center(),(real).1);
			iterate_cell(iter,fluid.mac_grid.grid){const VectorDi& cell=iter.Coord();
				VectorD pos=fluid.mac_grid.grid.Center(cell);
				fluid.levelset.phi(cell)=sphere.Phi(pos);}

		}break;}
		fluid.Enforce_Boundary_Conditions();}
};

#endif