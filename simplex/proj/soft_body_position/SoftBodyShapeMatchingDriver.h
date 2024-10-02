//#####################################################################
// Soft body FEM driver
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __SoftBodyShapeMatchingDriver_h__
#define __SoftBodyShapeMatchingDriver_h__
#include "Driver.h"
#include "SoftBodyShapeMatching.h"
#include "Constants.h"
#include "MeshFunc.h"

template<int d> class SoftBodyShapeMatchingDriver : public Driver
{Typedef_VectorDii(d);Typedef_MatrixD(d);using Base=Driver;
public:
	SoftBodyShapeMatching<d> soft_body;

	virtual void Initialize()
	{
		switch(test){
		case 1:{
			int m=15;int n=5;VectorD start=VectorD::Zero();real dx=(real)1/(real)m;
			soft_body.particles.Resize(m*n);
			int c=0;for(int i=0;i<m;i++)for(int j=0;j<n;j++){
				soft_body.particles.M(c)=(real)1;
				soft_body.particles.V(c)=VectorD::Zero();
				soft_body.particles.X(c++)=(start+AuxFunc::V<d>((real)i*dx,(real)j*dx));}
			MeshFunc::Rotate(soft_body.particles.XRef(),pi*(real).1);
		}break;}

		soft_body.Collision=std::bind(&SoftBodyShapeMatchingDriver::Collision,this,std::placeholders::_1);
		soft_body.Initialize();
	}

	void Kinematic_Boundary_Conditions(const real dt,const real time)
	{
	}

	void Collision(const real dt)
	{
		////impulse based
		//real ground=(real)-1.;real col_coef=(real)2e2;
		//for(int i=0;i<soft_body.particles.Size();i++){
		//	if(soft_body.particles.X(i)[1]<ground){
		//		soft_body.particles.F(i)+=col_coef*VectorD::Unit(1)*(ground-soft_body.particles.X(i)[1]);}}

		////position based
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

		{std::string file_name=frame_dir+"/particles";
		soft_body.particles.Write_To_File_3d(file_name);}

		std::cout<<"Write to frame "<<frame<<std::endl;
	}
};

#endif