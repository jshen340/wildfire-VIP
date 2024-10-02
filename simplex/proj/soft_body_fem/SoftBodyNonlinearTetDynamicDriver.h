//#####################################################################
// Soft body FEM driver
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __SoftBodyNonlinearTetDynamicDriver_h__
#define __SoftBodyNonlinearTetDynamicDriver_h__
#include <fstream>
#include "Common.h"
#include "Field.h"
#include "Driver.h"
#include "AuxFunc.h"
#include "MeshFunc.h"
#include "MeshFuncExt.h"
#include "SoftBodyNonlinearFemTet.h"

template<int d> class SoftBodyNonlinearTetDynamicDriver : public Driver
{Typedef_VectorDii(d);Typedef_MatrixD(d);using Base=Driver;
public:
	template<int d> using VolMesh=typename If<d==2,TriangleMesh<d>,TetrahedronMesh<d> >::Type;
	VolMesh<d> mesh;
	SoftBodyNonlinearFemTet<d> soft_body;

	virtual void Initialize()
	{
		soft_body.use_omp=true;
		//cfl=1e-5;
		cfl=1e-4;
		Initialize_Mesh();
		Initialize_Materials();
		Initialize_Boundary_Conditions();
	}

	void Initialize_Mesh()
	{
		VectorDi cell_counts=VectorDi::Ones();cell_counts*=scale;cell_counts[1]/=4;if(d==3)cell_counts[2]/=4;
		real dx=(real)1/cell_counts[0];
		switch(test){
		case 1:{
			MeshFunc::Initialize_Lattice_Mesh(cell_counts,dx,&mesh);	
		}break;
		case 2:{
			MeshFunc::Initialize_Lattice_Mesh(cell_counts,dx,&mesh);	
			soft_body.Kinematic_Boundary_Condition=std::bind(
				&SoftBodyNonlinearTetDynamicDriver::Kinematic_Boundary_Conditions,this,
				std::placeholders::_1,std::placeholders::_2);
			soft_body.use_body_force=false;
		}break;
		case 3:{
			if constexpr (d==2){
				MeshFunc::Initialize_Circle_Mesh(VectorD::Zero(),(real).5,&mesh,16);
				soft_body.Collision=std::bind(&SoftBodyNonlinearTetDynamicDriver::Collision,this,std::placeholders::_1);
				soft_body.use_damping=true;
				soft_body.damping*=(real)1;			
			}
			else std::cerr<<"Error: [SoftBodyNonlinearFemTet] dimension incompatable"<<std::endl;
		}break;
		case 4:{
			cell_counts[1]=cell_counts[0]/16;if(d==3)cell_counts[2]=cell_counts[1]/16;
			MeshFunc::Initialize_Lattice_Mesh(cell_counts,dx,&mesh);
		}break;
		}
		
		soft_body.Initialize(mesh);
	}

	void Initialize_Boundary_Conditions()
	{
		switch(test){
		case 1:{	////beam with one end fixed under gravity
			int axis=0;
			for(int i=0;i<(int)soft_body.mesh->Vertices().size();i++){
				const VectorD& pos=soft_body.mesh->Vertices()[i];
				if(pos[axis]==0){soft_body.Set_Fixed(i);}}
		}break;
		case 2:{	////beam with kinematic bc
			int axis=0;
			for(int i=0;i<(int)soft_body.mesh->Vertices().size();i++){
				const VectorD& pos=soft_body.mesh->Vertices()[i];
				if(pos[axis]==0){soft_body.Set_Displacement(i,VectorD::Unit(0)*(real)-1);}
				if(pos[axis]==1){soft_body.Set_Displacement(i,VectorD::Unit(0)*(real)1);}}
		}break;
		case 3:{	////a sphere falling onto the ground
			soft_body.use_body_force=true;
		}break;
		case 4:{	////multi-material
			int axis=0;
			for(int i=0;i<(int)soft_body.mesh->Vertices().size();i++){
				const VectorD& pos=soft_body.mesh->Vertices()[i];
				if(pos[axis]==0||pos[axis]==(real)1){soft_body.Set_Fixed(i);}}			
		}break;
		}
	}

	void Kinematic_Boundary_Conditions(const real dt,const real time)
	{
		switch(test){
		case 2:{
			real x=(real).1*sin(time*pi);
			for(auto& iter:soft_body.bc.psi_D_values){
				int i=iter.first;
				VectorD dis=iter.second;
				if(dis==VectorD::Unit(0)*(real)-1)soft_body.X()[i]=soft_body.X0[i]-VectorD::Unit(0)*x;
				else if(dis==VectorD::Unit(0)*(real)1)soft_body.X()[i]=soft_body.X0[i]+VectorD::Unit(0)*x;}
		}break;}
	}

	void Initialize_Materials()
	{
		soft_body.materials.clear();
		soft_body.materials.push_back(ElasticParam((real)1e2,(real).35));	

		AuxFunc::Fill(soft_body.material_id,0);
		switch(test){
		case 4:{	////multi-material
			soft_body.materials[0]=(ElasticParam((real)1,(real).35));
			soft_body.materials.push_back(ElasticParam((real)1e3,(real).35));	

			int s=(int)soft_body.material_id.size();
			int p=s/9;
			for(int i=0;i<s;i++){
				if((i/p)%2==0)soft_body.material_id[i]=0;
				else soft_body.material_id[i]=1;}
		}break;
		}
	}

	void Collision(const real dt)
	{
		real ground=(real)-1.;
		for(int i=0;i<soft_body.particles.Size();i++){
			if(soft_body.particles.X(i)[1]<ground){
				soft_body.particles.X(i)[1]=ground;
				soft_body.particles.V(i)[1]=(real)0;}}
	}

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		soft_body.Advance(dt,time);
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);

		{std::string file_name=frame_dir+(d==2?"/triangle_mesh":"/tetrahedron_mesh");
		soft_body.mesh->Write_To_File_3d(file_name);}

		{std::string file_name=frame_dir+"/particles";
		soft_body.particles.Write_To_File_3d(file_name);}

		{
			std::string file_name=frame_dir+"/mat";
			int n=(int)soft_body.material_id.size();
			Field<real,1> mat;mat.Resize(n);
			for(int i=0;i<n;i++){
				mat.array[i]=(real)soft_body.material_id[i];
				//std::cout<<mat.array[i]<<", ";
			}
			mat.Write_To_File_3d(file_name);
		}

		std::cout<<"Write to frame "<<frame<<std::endl;
	}
};

#endif