//#####################################################################
// Particle level set driver
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __ParticleLevelSetDriver_h__
#define __ParticleLevelSetDriver_h__
#include "Common.h"
#include "File.h"
#include "ParticleLevelSet.h"
#include "GeometryPrimitives.h"
#include "AnalyticalFields.h"
#include "MarchingCubes.h"
#include "Driver.h"

template<int d> class ParticleLevelSetDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;
public:
	real domain_size=(real)1;
	std::shared_ptr<VorticityField<d> > velocity_field=nullptr;
	ParticleLevelSet<d> pls;

	virtual void Initialize()
	{
		VectorDi cell_counts=VectorDi::Ones()*scale;real dx=domain_size/(real)scale;
		VectorD domain_min=VectorD::Zero();
		Grid<d> grid(cell_counts,dx,domain_min);
		pls.Initialize(grid);
		
		Sphere<d> sphere(grid.Center()+VectorD::Unit(1)*domain_size*(real).25,domain_size*(real).1);
		iterate_cell(iter,pls.grid){const VectorDi& cell=iter.Coord();
			pls.phi(cell)=sphere.Phi(pls.grid.Center(cell));}
		pls.Initialize_Particles();
		pls.Fast_Marching();

		velocity_field=std::make_shared<VorticityField<d> >();
		velocity_field->origin=grid.Center();
		velocity_field->test=1;	
	}
	
	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		Field<real,d> ghost_phi=pls.phi;
		Interpolation<d> intp(pls.grid);
		////advect levelset
		iterate_cell(iter,pls.grid){const VectorDi& cell=iter.Coord();
			VectorD p1=pls.grid.Center(cell);
			VectorD v=velocity_field->Velocity(p1);
			VectorD p2=pls.grid.Center(cell)-v*dt*(real).5;
			v=velocity_field->Velocity(p2);
			VectorD x=p1-v*dt;
			pls.phi(cell)=intp.Interpolate_Centers(ghost_phi,x);}
		////advect particles
		for(int i=0;i<pls.particles.Size();i++){
			VectorD v=velocity_field->Velocity(pls.particles.X(i));
			VectorD p2=pls.particles.X(i)+v*dt*(real).5;
			v=velocity_field->Velocity(p2);
			pls.particles.X(i)+=v*dt;}
		////particle correction and fast marching
		//pls.Correction();
		pls.Fast_Marching();
		//pls.Correction();
		//pls.Update_Particle_Radii();
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);
		if(frame==0){
			std::string file_name=frame_dir+"/grid";
			pls.grid.Write_To_File_3d(file_name);
			std::cout<<"Write to file "<<file_name<<std::endl;}
		
		{std::string file_name=frame_dir+"/phi";
		pls.phi.Write_To_File_3d(file_name);}

		{std::string file_name=frame_dir+"/particles";
		pls.particles.Write_To_File_3d(file_name);}

		MarchingCubes<d> marching_cubes(pls);
		marching_cubes.Marching();
		if(d==2){
			std::string file_name=frame_dir+"/segment_mesh";
			(*marching_cubes.mesh).Write_To_File_3d(file_name);}
		else{
			std::string file_name=frame_dir+"/triangle_mesh";
			(*marching_cubes.mesh).Write_To_File_3d(file_name);}

		{std::string file_name=frame_dir+"/curvature";
		Field<real,d> curv(pls.grid.cell_counts,(real)0);
		iterate_cell(iter,pls.grid){const VectorDi& cell=iter.Coord();
			if(abs(pls.phi(cell))>(real)pls.grid.dx*(real)5)continue;
			curv(cell)=pls.Curvature(pls.grid.Center(cell));}
		curv.Write_To_File_3d(file_name);}

		{std::string file_name=frame_dir+"/normal";
		Field<VectorD,d> normal(pls.grid.cell_counts,VectorD::Zero());
		iterate_cell(iter,pls.grid){const VectorDi& cell=iter.Coord();
			if(abs(pls.phi(cell))>(real)pls.grid.dx*(real)5)continue;
			normal(cell)=pls.Normal(pls.grid.Center(cell));}
		Field<Vector3,3> n3;VF_Dim_Conversion<real,d,3>(normal,n3);
		n3.Write_To_File_3d(file_name);}

		std::cout<<"Write to frame "<<frame<<std::endl;
	}
};

#endif