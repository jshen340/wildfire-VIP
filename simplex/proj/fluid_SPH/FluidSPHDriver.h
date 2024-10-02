//#####################################################################
// SPH Fluid driver
// Copyright (c) (2018-), Xiangxin Kong
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//#####################################################################

#ifndef __FluidSPHDriver_h__
#define __FluidSPHDriver_h__

#include "Common.h"
#include "File.h"
#include "GeometryPrimitives.h"
#include "Driver.h"
#include "FluidSPH.h"
#include "NeighborSearcher.h"
#include "RenderFunc.h"

template<int d> class FluidSPHDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;
public:
	FluidSPHND<d> fluid;

	virtual real CFL() const
	{
		return cfl*fluid.CFL();
	}

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		fluid.Advance(dt);
	}

	virtual void Write_Output_Files(const int frame)
	{
		Base::Write_Output_Files(frame);

		////SPH particle X, F, V
		RenderFunc::Write_Points_Float<d, real>(frame_dir + "/tracker_points", fluid.particles.XRef());
		RenderFunc::Write_Vectors_Float<d, real>(frame_dir + "/point_velocity", fluid.particles.XRef(), fluid.particles.VRef());
		RenderFunc::Write_Vectors_Float<d, real>(frame_dir + "/point_force", fluid.particles.XRef(), fluid.particles.FRef());			
		std::cout<<"write frame "<<frame<<std::endl;
	}

	virtual void Initialize()
	{
		switch (test) {
		case 1: Case_1(); break; ////droplet with central gravity
		case 2: Case_2(); break; ////dam break
		case 3: Case_3(); break; //test for packing algorithm
		}
	
	}

	void Collision(const real dt)
	{
		switch(test){
		case 2:{	////tank
			int pn=fluid.particles.Size();
			for(int i=0;i<pn;i++){
				if(fluid.particles.X(i)[0]>(real)1){
					fluid.particles.X(i)[0]=(real)1;
					fluid.particles.V(i)[0]=(real)0;}
				else if(fluid.particles.X(i)[0]<(real)0){
					fluid.particles.X(i)[0]=(real)0;
					fluid.particles.V(i)[0]=(real)0;}
				if(fluid.particles.X(i)[1]<(real)0){
					fluid.particles.X(i)[1]=(real)0;
					fluid.particles.V(i)[1]=(real)0;}}		
		}break;
		}
	}

	void Case_1(void);
	void Case_2(void);
	void Case_3(void);

protected:
	void Initialize_Particle(const int i,const VectorD& pos)
	{
		fluid.particles.X(i)=pos;
		fluid.particles.V(i)=VectorD::Zero();
		fluid.particles.F(i)=VectorD::Zero();
		fluid.particles.M(i)=(real)1;
		fluid.particles.C(i)=1.;
		fluid.particles.I(i)=i;
	}
};

#endif