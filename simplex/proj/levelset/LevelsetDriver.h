//////////////////////////////////////////////////////////////////////////
// Level set driver
// Bo Zhu
//////////////////////////////////////////////////////////////////////////
#ifndef __LevelSetDriver_h__
#define __LevelSetDriver_h__
#include "Common.h"
#include "File.h"
#include "LevelSet.h"
#include "GeometryPrimitives.h"
#include "MarchingCubes.h"

template<int d> class LevelSetDriver
{Typedef_VectorDii(d);
public:
	int test=1;
	std::string output_directory="output";
	int scale=8;
	real domain_size=(real)1;
	LevelSet<d> levelset;

	LevelSetDriver(){}

	void Initialize()
	{
		VectorDi cell_counts=VectorDi::Ones()*scale;real dx=domain_size/(real)scale;
		VectorD domain_min=-VectorD::Ones()*domain_size*(real).5;
		Grid<d> grid(cell_counts,dx,domain_min);levelset.Initialize(grid);
		Sphere<d> sphere(VectorD::Zero(),domain_size*(real).2);
		iterate_cell(iter,levelset.grid){const VectorDi& cell=iter.Coord();
			levelset.phi(cell)=sphere.Phi(levelset.grid.Center(cell));}
		levelset.Fast_Marching();
	}
	
	void Run()
	{
		Write_Output_Files(0);
	}

	void Write_Output_Files(const int frame)
	{	
		if(!File::Directory_Exists(output_directory.c_str()))File::Create_Directory(output_directory);
		std::string frame_dir=output_directory+"/"+std::to_string(frame);
		if(!File::Directory_Exists(frame_dir.c_str()))File::Create_Directory(frame_dir);

		{std::string file_name=frame_dir+"/grid";
		levelset.grid.Write_To_File_3d(file_name);}
		
		{std::string file_name=frame_dir+"/phi";
		levelset.phi.Write_To_File_3d(file_name);}

		MarchingCubes<d> marching_cubes(levelset);
		marching_cubes.Marching();
		if(d==2){
			std::string file_name=frame_dir+"/segment_mesh";
			(*marching_cubes.mesh).Write_To_File_3d(file_name);}
		else{
			std::string file_name=frame_dir+"/triangle_mesh";
			(*marching_cubes.mesh).Write_To_File_3d(file_name);}

		{std::string file_name=frame_dir+"/curvature";
		Field<real,d> curv(levelset.grid.cell_counts,(real)0);
		iterate_cell(iter,levelset.grid){const VectorDi& cell=iter.Coord();
			if(abs(levelset.phi(cell))>(real)levelset.grid.dx*(real)5)continue;
			curv(cell)=levelset.Curvature(levelset.grid.Center(cell));
			std::cout<<curv(cell)<<", ";}
		curv.Write_To_File_3d(file_name);}

		{std::string file_name=frame_dir+"/normal";
		Field<VectorD,d> normal(levelset.grid.cell_counts,VectorD::Zero());
		iterate_cell(iter,levelset.grid){const VectorDi& cell=iter.Coord();
			if(abs(levelset.phi(cell))>(real)levelset.grid.dx*(real)5)continue;
			normal(cell)=levelset.Normal(levelset.grid.Center(cell));}
		Field<Vector3,3> n3;VF_Dim_Conversion<real,d,3>(normal,n3);
		n3.Write_To_File_3d(file_name);}
	}
};

#endif