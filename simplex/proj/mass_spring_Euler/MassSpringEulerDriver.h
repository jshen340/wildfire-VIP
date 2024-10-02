//#####################################################################
// Mass spring driver
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __MassSpringEulerDriver_h__
#define __MassSpringEulerDriver_h__
#include <memory>
#include "Common.h"
#include "Driver.h"
#include "MassSpringEuler.h"
#include "RandomNumber.h"

template<int d> class MassSpringEulerDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;using Base::verbose;
public:
	MassSpringEuler<d> system;
	//int times;

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		system.Advance(dt);
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);
		
		
		////Write grid
		{
			std::string file_name = frame_dir + "/grid";
			system.grid.Write_To_File_3d(file_name);
		}

		////Write potential
		{
			std::string file_name = frame_dir + "/phi";
			system.potential.Write_To_File_3d(file_name);
		}

		////Write Particles
		{std::string file_name=frame_dir+"/particles";
		system.particles.Write_To_File_3d(file_name);}

		////Write BC
		{std::string file_name=frame_dir+"/psi_D";
		Particles<d> psi_D_particles;
		for(auto p:system.psi_D_values){int idx=p.first;
			int i=psi_D_particles.Add_Element();psi_D_particles.X(i)=system.particles.X(idx);}
		psi_D_particles.Write_To_File_3d(file_name);}
		
		
		

		std::ofstream out_file;
		std::string file_name = frame_dir + "/../" + std::to_string(frame) + ".csv";
		out_file.open(file_name);


		for (int i = 0; i < system.particles.Size(); i++) {
			//std::cout << i << planets.particles.X(i) << ", Vi, " << planets.particles.V(i) << ", " << std::endl;
			out_file << system.particles.M(i) << std::endl;
			out_file << system.particles.X(i) << std::endl << system.particles.V(i) << std::endl;
		}

		out_file.close();
	}
	
	virtual void Initialize()
	{
		verbose=true;

		switch(test){
		case 1:{
			system.use_body_force = false;
			int n = 100;
			system.particles.Resize(n);
			VectorDi cell_counts = VectorDi::Ones() * 128; real dx = (real)2 / ((real)cell_counts[0]);
			system.grid.Initialize(cell_counts, dx, -VectorD::Ones());

			real offset = (system.grid.domain_max[0] - system.grid.domain_min[0]) * (real).2;
			RandomNumber rand(system.grid.domain_min[0] + offset, system.grid.domain_max[0] - offset);// RandomNumber does not support times
			for(int i=0;i<n;i++){
				VectorD x = rand.VectorValue<d>();
				system.particles.X(i) = x;
				system.particles.M(i)=(real)1;
				//std::cout << x.transpose() << std::endl;
			}
			//std::cout << "next" << std::endl;
		}break;
		}

		system.Initialize();
	}
};

#endif