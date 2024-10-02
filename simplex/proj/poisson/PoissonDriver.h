//#####################################################################
// Poisson driver
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __PoissonDriver_h__
#define __PoissonDriver_h__
#include "Common.h"
#include "File.h"
#include "PoissonOld.h"
#include "PoissonMtxFreeOld.h"

template<int d> class PoissonDriver
{Typedef_VectorDii(d);
public:
	int test=1;
	std::string output_directory="output";
	int scale=8;
	real domain_size=(real)1;
	PoissonOld<d> poisson;
	//PoissonMtxFreeOld<d> poisson;

	PoissonDriver(){}

	void Initialize()
	{
		VectorDi cell_counts=VectorDi::Ones()*scale;real dx=domain_size/(real)scale;
		VectorD domain_min=-VectorD::Ones()*domain_size*(real).5;
		poisson.Initialize(cell_counts,dx,domain_min);
	}
	
	void Run()
	{
		Update_Source();
		Update_Boundary_Conditions();
		poisson.Build_And_Solve();
		Write_Output_Files(0);
	}

	void Update_Source()
	{
		iterate_cell(iter,poisson.mac_grid.grid){
			const VectorDi& cell=iter.Coord();
			const VectorD center=poisson.mac_grid.grid.Center(cell);
			poisson.rhs(cell)=lap_F(center,test);}
	}

	void Update_Boundary_Conditions()
	{
		if (test==4) {
			//set Dirichlet boundary on the other sides
			iterate_cell(iter,poisson.mac_grid.grid){const VectorDi& cell=iter.Coord();
			if(poisson.mac_grid.grid.Is_Boundary_Cell(cell)&&cell[1]!=0){
				const VectorD center=poisson.mac_grid.grid.Center(cell);
				real psi_D_value=(real)0;
			}}

			//set Neuman boundary on the bottom border
			for (int axis = 0; axis < d; axis++) {
				int face_num = poisson.mac_grid.face_grids[axis].Number_Of_Nodes(); 
				for (int i = 0; i < face_num; i++) {
					VectorDi face = poisson.mac_grid.face_grids[axis].Node_Coord(i);
					VectorD pos = poisson.mac_grid.Face_Center(axis, face);
					if (poisson.mac_grid.Is_Axial_Boundary_Face(axis, face)&&face[1]==0) {
						poisson.bc.Set_Psi_N(axis,face,(real)1);
					}
				}
			}
		}
		else {
			iterate_cell(iter,poisson.mac_grid.grid){const VectorDi& cell=iter.Coord();
			if(poisson.mac_grid.grid.Is_Boundary_Cell(cell)){
				const VectorD center=poisson.mac_grid.grid.Center(cell);
				real psi_D_value=F(center,test);
				//real psi_D_value=(real)0;
				//poisson.Set_Psi_D(cell,psi_D_value);	////dbg
				//poisson.bc.Set_Psi_D(cell,(ushort)CellType::Fluid,psi_D_value);
			}}
		}
	}

	void Write_Output_Files(const int frame)
	{	
		if(!File::Directory_Exists(output_directory.c_str()))File::Create_Directory(output_directory);
		std::string frame_dir=output_directory+"/"+std::to_string(frame);
		if(!File::Directory_Exists(frame_dir.c_str()))File::Create_Directory(frame_dir);

		{std::string file_name=frame_dir+"/grid";
		poisson.mac_grid.grid.Write_To_File_3d(file_name);}
		
		{std::string file_name=frame_dir+"/x";
		poisson.p.Write_To_File_3d(file_name);}
	}

	////test functions
	////test 1: f=sin(x)
	////test 2: f=x^2+y^2
	////test 3: f=x^2+y^2+z^2
	real F(const VectorD& p,const int test)
	{
		switch(test){
		case 1:return sin(p[0]);
		case 2:return pow(p[0],2)+pow(p[1],2);
		case 3:return pow(p[0],2)+pow(p[1],2)+pow(p[2],2);
		default:return (real)0;break;}
	}

	VectorD dF(const VectorD& p)
	{
		switch(test){ 
		case 1:return AuxFunc::V<d>(cos(p[0]), (real)0, (real)0);
		case 2:return AuxFunc::V<d>(2*p[0], 2*p[1], (real)0);
		case 3:return AuxFunc::V<d>(2*p[0], 2*p[1], 2*p[2]);
		default:return VectorD::Zero();break;}
	}

	real lap_F(const VectorD& p,const int test)
	{
		switch(test){
		case 1:return -sin(p[0]);
		case 2:return (real)4;
		case 3:return (real)6;
		default:return (real)0;}
	}
};

#endif