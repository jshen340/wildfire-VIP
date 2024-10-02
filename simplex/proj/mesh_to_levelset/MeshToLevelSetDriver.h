//////////////////////////////////////////////////////////////////////////
// Mesh To Level Set Driver
// Copyright (c) (2020-), Fan Feng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __MeshToLevelSetDriver_h__
#define __MeshToLevelSetDriver_h__
#include "Common.h"
#include "Driver.h"
#include "Field.h"
#include "Grid.h"
#include "File.h"
#include "Mesh.h"
#include "MeshToLevelSet.h"
#include "MeshAdvFunc.h"
#include "MarchingCubes.h"
#include <iostream>
#include <fstream>
#ifdef USE_TINY_OBJ_LOADER
#include "TinyObjLoader.h"
#endif

template<int d,bool use_mat_id=false> class MeshToLevelSetDriver : public Driver
{Typedef_VectorDii(d);Typedef_MatrixD(d);using Base=Driver;
public:
	Grid<d> grid;
	Field<real,d> rho;
	Field<real,d> rho2;
	SurfaceMesh<d> mesh;
	LevelSet<d> level_set;
	std::string input_path;

	void Initialize()
	{
		VectorDi cell_counts=VectorDi::Ones()*scale;
		real dx=(real)1/(real)scale;
		grid.Initialize(cell_counts,dx);
		rho.Resize(grid.cell_counts,(real)0);
		rho2.Resize(grid.cell_counts,(real)0);

		if (!input_path.empty()){std::cout << "File path is " << input_path << std::endl;}
		else{std::cerr<<"Please specify input path!"<<std::endl; }
		
		switch(test){
			case 1:{
			if constexpr (d == 2) {
				////read from mesh file
				std::shared_ptr<VolumetricMesh<2> > volume_mesh = std::make_shared<VolumetricMesh<2> >();
				MeshFunc::Initialize_Polygon_Mesh_From_File(input_path.c_str(), volume_mesh.get());
				mesh.vertices = volume_mesh->vertices;
				MeshFunc::Volumetric_Mesh_Surface(volume_mesh->Elements(),mesh.Elements());
			}else if constexpr (d == 3) {
				////load 3d obj mesh and convert to txt
				/*Array<std::shared_ptr<TriangleMesh<3> > > tri_meshes;
				Obj::Read_From_Obj_File(input_path, tri_meshes);
				mesh = *tri_meshes[0];
				std::string output_path="C:\\Users\\fan\\Research\\simplex\\data\\meshes\\rabbit_scaled_complete.txt";
				std::ofstream output(output_path);
				File::Write_Text(output,mesh);*/
				////read in txt file
				std::ifstream input(input_path.c_str());
				File::Read_Text(input,mesh);
			}}break;
		}
	}

	void Run()
	{
		level_set.Initialize(grid);
		MeshFunc::Mesh_To_Level_Set(mesh,grid,level_set);
		std::cout<<"converted mesh to levelset successfully!"<<std::endl;

		int cell_num=grid.cell_counts.prod();
		for(int i=0;i<cell_num;i++){
			const VectorDi cell=grid.Cell_Coord(i);
			rho(cell)=(real)level_set.phi(cell);}

		for(int i=0;i<cell_num;i++){
			const VectorDi cell=grid.Cell_Coord(i);
			rho2(cell)=(real)level_set.phi(cell);}
		Write_Output_Files(0);
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);

		if(frame==0){
			std::string file_name=frame_dir+"/grid";
			grid.Write_To_File_3d(file_name);}

		{std::string file_name=frame_dir+"/rho";
		rho.Write_To_File_3d(file_name);}

		{std::string file_name=frame_dir+"/rho2";
		rho2.Write_To_File_3d(file_name);}

		{MarchingCubes<d> marching_cubes(level_set);
		marching_cubes.Marching();
		if(d==2){
			std::string file_name=frame_dir+"/segment_mesh";
			if((*marching_cubes.mesh).Vertices().size()>0)
				(*marching_cubes.mesh).Write_To_File_3d(file_name);}
		else{
			std::string file_name=frame_dir+"/triangle_mesh";
			if((*marching_cubes.mesh).Vertices().size()>0)
				(*marching_cubes.mesh).Write_To_File_3d(file_name);}}

		{std::string file_name=frame_dir+"/phi_"+std::to_string(scale);
		File::Write_Binary_To_File(file_name,level_set.phi);}
	}
};
#endif