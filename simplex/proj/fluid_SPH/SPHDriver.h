//#####################################################################
// Universal SPH Driver
// Copyright (c) (2018-), Xiangxin Kong, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//#####################################################################

#ifndef __SPHDriver_h__
#define __SPHDriver_h__

#include<iostream>

#include "Common.h"
#include "File.h"
#include "GeometryPrimitives.h"
#include "FluidSPH.h"
#include "Particles.h"
#include "Driver.h"
#include "RenderFunc.h"

template<int d> class SPHDriver : public Driver
{Typedef_VectorDii(d); using Base = Driver;
public:
	//0:weakly compressible, 1:PCISPH, 2:IISPH
	int algo_type = 0;
	std::shared_ptr<FluidSPH<d> > fluid;

	virtual real CFL(void) const {
		real cfl_time = fluid->dx / fluid->Max_Velocity() * Base::cfl;
		return AuxFunc::Clamp(cfl_time, 0.01 / Base::frame_rate, 1.0 / Base::frame_rate);
	}

	virtual void Advance_One_Time_Step(const real dt, const real time)
	{
		Timer timer;
		fluid->Advance(dt, time);
		int iter_per_frame = (int)(1.0 / Base::frame_rate / dt + 0.5);
		if (dt < 1e-8) iter_per_frame = 0;
		if (verbose) std::cout << "["
			<< std::setw(6) << iter_per_frame << " Iterations Per Frame]     ... "
			<< std::setw(7) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << timer.Elapse(PhysicalUnits::s) << "s used, current time: " << time << std::endl;
		std::cout << std::defaultfloat;
	}

	virtual void Write_Output_Files(const int frame)
	{
		Base::Write_Output_Files(frame);
		////SPH particle X, F, V
		RenderFunc::Write_Points_Float<d, real>(frame_dir + "/tracker_points", fluid->particles.XRef());
		RenderFunc::Write_Vectors_Float<d, real>(frame_dir + "/point_velocity", fluid->particles.XRef(), fluid->particles.VRef());
	}

	virtual void Initialize()
	{
		Base::verbose = true;
		cfl = 0.5;
		switch (test) {
		case 1: Case_1(); break;//a sphere drop inside a sphere bowl
		case 2: Case_2(); break;//a cubic drop inside a cubic bounding bow
		case 3: Case_3(); break;//dam break
		case 4: Case_4(); break;//hydrostatic, water fills the lower part of cubic
		}
		std::cout << "#     Initialized " << fluid->particles.Size() << " particles\n";
	}

	void Add_SPH_Particle(SPHParticles<d>& particles, const VectorD& pos, const real& m);
	void Case_1(void);
	void Case_2(void);
	void Case_3(void);
	void Case_4(void);

	protected:
};

#endif

