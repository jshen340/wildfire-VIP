//#####################################################################
// Poisson driver
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __PoissonIrregularDriver_h__
#define __PoissonIrregularDriver_h__
#include "Common.h"
#include "File.h"
#include "PoissonIrregularOld.h"
#include "MarchingCubes.h"
#include "Levelset.h"
#include "AuxFunc.h"

template<int d> class PoissonIrregularDriver
{Typedef_VectorDii(d);
public:
	int test=1;
	int boundary=1;
	int boundary_type=1;
	std::string output_directory="output";
	int scale=8;
	real domain_size=(real)1;
	PoissonIrregularOld<d> poisson;

	PoissonIrregularDriver(){}

	void Initialize()
	{
		VectorDi cell_counts=VectorDi::Ones()*scale;real dx=domain_size/(real)scale;
		VectorD domain_min=-VectorD::Ones()*domain_size*(real).5;
		MacGrid<d> mac_grid;mac_grid.Initialize(cell_counts,dx,domain_min);
		poisson.Initialize(mac_grid);
		std::cout<<"eps is: "<<poisson.eps<<std::endl;

		switch(boundary){
		case 1:{ //half irregular and half regular boundary
			VectorD center;real r;
			center=mac_grid.grid.Center();center[1]=domain_min[1];r=domain_size*(real).6;
			Sphere<d> sphere(center,r);
			iterate_cell(iter,mac_grid.grid){const VectorDi& cell=iter.Coord();
			poisson.levelset->phi(cell)=sphere.Phi(mac_grid.grid.Center(cell));}
		}break;
		case 2:{ //all irregular boundary
			VectorD center;real r;
			center=mac_grid.grid.Center();r=domain_size*(real).3;
			Sphere<d> sphere(center,r);
			iterate_cell(iter,mac_grid.grid){const VectorDi& cell=iter.Coord();
			poisson.levelset->phi(cell)=sphere.Phi(mac_grid.grid.Center(cell));}
		}break;
		case 3: { //a submerged circle and an halfly submerged sphere
			VectorD center1, center2; real r1, r2;
			center1=mac_grid.grid.Center();center1[1]+=(center1[1]-domain_min[1])*0.5;r1=domain_size*(real).2;
			center2=mac_grid.grid.Center();center2[1]=domain_min[1];r2=domain_size*(real).5;
			Sphere<d> sphere1(center1,r1); Sphere<d> sphere2(center2,r2);
			iterate_cell(iter,mac_grid.grid){const VectorDi& cell=iter.Coord();
			poisson.levelset->phi(cell)=std::min(sphere1.Phi(mac_grid.grid.Center(cell)), sphere2.Phi(mac_grid.grid.Center(cell)));}
		}break;
		}	
		
		poisson.psi_D_irregular = Irregular_Boundary_Psi_D();
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
		poisson.rhs.Fill((real)0);
		iterate_cell(iter,poisson.mac_grid.grid){
			const VectorDi& cell=iter.Coord();
			const VectorD pos=poisson.mac_grid.grid.Center(cell);
			if(poisson.Is_Inside_Levelset_eps(cell))
				poisson.rhs(cell)=lap_F(pos);}
	}

	void Update_Boundary_Conditions()
	{
		if (boundary_type == 1) {
			for (int axis = 0; axis < d; axis++) {
				int face_num = poisson.mac_grid.face_grids[axis].Number_Of_Nodes(); 
				for (int i = 0; i < face_num; i++) {
					VectorDi face = poisson.mac_grid.face_grids[axis].Node_Coord(i);
					VectorD pos = poisson.mac_grid.Face_Center(axis, face);
					if (poisson.mac_grid.Is_Axial_Boundary_Face(axis, face) && poisson.Is_Inside_Levelset_eps(pos)) {
						poisson.bc->Set_Psi_N(axis,face,dF(pos)[axis]);
					}
				}
			}
		}
		else if (boundary_type == 2) {
			iterate_cell(iter,poisson.mac_grid.grid){const VectorDi& cell=iter.Coord();
				if(poisson.mac_grid.grid.Is_Boundary_Cell(cell) && poisson.Is_Inside_Levelset_eps(cell)){
					const VectorD center=poisson.mac_grid.grid.Center(cell);
					real psi_D_value=F(center);poisson.bc->Set_Psi_D(cell,(ushort)CellType::Fluid,psi_D_value);
					real phi=poisson.levelset->Phi(center);
					real phi2=poisson.levelset->phi(cell);}}
		}
		else {
			std::cout << "Invalid boundary type: " << boundary_type << std::endl;
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

		////Write interface
		std::string file_name=frame_dir+"/phi";
		poisson.levelset->phi.Write_To_File_3d(file_name);

		MarchingCubes<d> marching_cubes((*poisson.levelset));
		marching_cubes.Marching();
		if(d==2){
			std::string file_name=frame_dir+"/segment_mesh";
			(*marching_cubes.mesh).Write_To_File_3d(file_name);}
		else{
			std::string file_name=frame_dir+"/triangle_mesh";
			(*marching_cubes.mesh).Write_To_File_3d(file_name);}
	}

	////test functions
	////test 1: f=sin(x)
	////test 2: f=x^2+y^2
	////test 3: f=x^2+y^2+z^2
	real F(const VectorD& p)
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

	real lap_F(const VectorD& p)
	{
		switch(test){
		case 1:return -sin(p[0]);
		case 2:return (real)4;
		case 3:return (real)6;
		default:return (real)0;}
	}

	std::function<real(const VectorD&)> Irregular_Boundary_Psi_D()
	{
		switch (test) {
		case 1: return [=](const VectorD& p)->real{return sin(p[0]);}; 
		case 2: return [=](const VectorD& p)->real{return pow(p[0],2)+pow(p[1],2);}; 
		case 3: return [=](const VectorD& p)->real{return pow(p[0],2)+pow(p[1],2)+pow(p[2],2);};
		default: [=](const VectorD& p)->real{return (real)FLT_MAX;};}
	}
};

#endif