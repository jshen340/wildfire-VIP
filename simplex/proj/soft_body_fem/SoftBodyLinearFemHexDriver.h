//#####################################################################
// Soft body FEM driver
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __SoftBodyLinearFemHexDriver_h__
#define __SoftBodyLinearFemHexDriver_h__
#include <fstream>
#include "Common.h"
#include "Driver.h"
#include "SoftBodyLinearFemHex.h"

template<int d> class SoftBodyLinearFemHexDriver : public Driver
{Typedef_VectorDii(d);Typedef_MatrixD(d);using Base=Driver;
public:
	HexMesh<d> mesh;
	SoftBodyLinearFemHex<d> soft_body;

	void Initialize()
	{
		Initialize_Mesh();
		Initialize_Materials();
		Initialize_Boundary_Conditions();
	}

	void Initialize_Mesh()
	{
		VectorDi cell_counts=VectorDi::Ones();cell_counts*=scale;//cell_counts[1]/=2;
		real dx=(real)1/cell_counts[0];
		switch(test){
		case 1:{
			MeshFunc::Initialize_Lattice_Mesh<d>(cell_counts,dx,&mesh);			
		}break;
		case 2:{
			MeshFunc::Initialize_Lattice_Mesh<d>(cell_counts,dx,&mesh);	
		}break;}
		
		soft_body.Initialize(mesh);
	}

	void Initialize_Materials()
	{
		soft_body.materials.clear();
		soft_body.materials.push_back(ElasticParam((real)1,(real).3));	
		AuxFunc::Fill(soft_body.material_id,0);
	}

	void Initialize_Boundary_Conditions()
	{
		switch(test){
		case 1:{
			int axis=0;VectorD dis=VectorD::Zero();dis[axis]=(real)1;
			for(int i=0;i<(int)soft_body.mesh->Vertices().size();i++){
				const VectorD& pos=soft_body.mesh->Vertices()[i];
				if(pos[axis]==0){soft_body.Set_Fixed(i);std::cout<<"fix "<<i<<": "<<pos.transpose()<<std::endl;}
				else if(pos[axis]==1){soft_body.Set_Displacement(i,dis);std::cout<<"dis "<<i<<": "<<pos.transpose()<<std::endl;}}
		}break;
		case 2:{
			VectorD f=VectorD::Zero();f[1]=(real)-1;
			for(int i=0;i<(int)soft_body.mesh->Vertices().size();i++){
				const VectorD& pos=soft_body.mesh->Vertices()[i];
				if(pos[0]==0){soft_body.Set_Fixed(i);}
				else if(pos[0]==1&&pos[1]==0){soft_body.Add_Force(i,f);}}
		}break;}
	}

	void Run()
	{
		soft_body.Allocate_K();
		soft_body.Update_K_And_f();
		soft_body.Solve();

		Write_Output_Files(0);
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);

		if(frame==0){
			std::string file_name=frame_dir+(d==2?"/quad_mesh":"/hex_mesh");
			soft_body.mesh->Write_To_File_3d(file_name);
			std::cout<<"Write to file "<<file_name<<std::endl;}
		
		{std::string file_name=frame_dir+"/displacement";
		Write_To_Field1_3d<VectorX,real,d>(soft_body.u,file_name);}

		{std::string file_name=frame_dir+"/force";
		Write_To_Field1_3d<VectorX,real,d>(soft_body.f,file_name);}

		//{std::string file_name=frame_dir+"/stress";
		//Field<MatrixD,1> stress;soft_body.Compute_Stress(stress.array);
		//file_name=frame_dir+"/strain";
		//Field<MatrixD,1> strain;soft_body.Compute_Strain(strain.array);

		//Field<real,1> von_mises;von_mises.Resize((int)stress.array.size());
		//for(auto i=0;i<(int)stress.array.size();i++){von_mises.array[i]=LinearFemFunc<d>::Von_Mises_Stress(stress.array[i]);}
		//file_name=frame_dir+"/von_mises_stress";}
	}
};

#endif