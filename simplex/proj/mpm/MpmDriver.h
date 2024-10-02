//#####################################################################
// Material point method driver
// Copyright (c) (2018-), Jiayin Hu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __MpmDriver_h__
#define __MpmDriver_h__
#include "Common.h"
#include "File.h"
#include "GeometryPrimitives.h"
#include "MaterialPoint.h"
#include "Particles.h"
#include "Driver.h"
#include "RandomNumber.h"

template<int d> class MpmDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;
public:
	MaterialPoint<d> material;

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		material.Advance(dt);
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);

		if(frame==0){
			std::string file_name = frame_dir + "/grid";
			material.grid.Write_To_File_3d(file_name);}
		{std::string file_name = frame_dir + "/particles";
		material.particles.Write_To_File_3d(file_name);}
		{std::string file_name = frame_dir + "/velocity";
		material.velocity.Write_To_File_3d(file_name);}
		std::cout<<"Write files to frame "<<frame<<std::endl;
	}
	
	virtual void Initialize()
	{
		real length = (real)1;
		VectorDi cell_counts = VectorDi::Ones() * 80;
		real dx = length / (real)cell_counts[0];
		switch(test){
		case 1:{
			material.use_body_force = true;
			int n = 500;
			material.particles.Resize(3*n);
			material.grid.Initialize(cell_counts, dx, VectorD::Zero());

			RandomNumber rand(-0.08, 0.08);
			VectorD origin = VectorD::Zero();
			origin[0] = 0.55; origin[1] = 0.45;
			VectorD origin2 = VectorD::Zero();
			origin2[0] = 0.45; origin2[1] = 0.65;
			VectorD origin3 = VectorD::Zero();
			origin3[0] = 0.55; origin3[1] = 0.85;
			for (int i = 0; i < n; i++) {
				////object-0
				material.particles.M(i) = (real)1;
				material.particles.V(i) = VectorD::Zero();
				material.particles.X(i) = origin + rand.VectorValue<d>();
				////object-1
				material.particles.M(n+i) = (real)1;
				material.particles.V(n+i) = VectorD::Zero();
				material.particles.X(n+i) = origin2 + rand.VectorValue<d>();
				////object-2
				material.particles.M(2*n + i) = (real)1;
				material.particles.V(2*n + i) = VectorD::Zero();
				material.particles.X(2*n + i) = origin3 + rand.VectorValue<d>();}
		}break;}
		material.Initialize();
	}
};

#endif