//#####################################################################
// Levelset main
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __Main_cpp__
#define __Main_cpp__
#include <iostream>
#include "ParseArgs.h"
#include "LevelsetDriver.h"
#include "ParticleLevelsetDriver.h"
#include "ParticleLevelsetAdvectionDriver.h"
#include "FluidParticleLevelsetDriver.h"

int main(int argc,char* argv[])
{
	const int d=2;
    ParseArgs parse_args;
    parse_args.Add_String_Argument("-o","output","output path");
    parse_args.Add_Integer_Argument("-s",64,"resolution");
	parse_args.Add_Integer_Argument("-test",1,"test");
	parse_args.Add_Integer_Argument("-driver",1,"driver");
	parse_args.Add_Integer_Argument("-use_pls",1,"use particle levelset?1:0");

    parse_args.Parse(argc,argv);
	int driver_id=parse_args.Get_Integer_Value("-driver");
	switch(driver_id){
	case 1:{
		LevelSetDriver<d> driver;
		driver.output_directory=parse_args.Get_String_Value("-o");
		driver.scale=parse_args.Get_Integer_Value("-s");
		driver.test=parse_args.Get_Integer_Value("-test");

		driver.Initialize();
		driver.Run();	
	}break;
	case 2:{
		ParticleLevelSetDriver<d> driver;
		driver.output_dir=parse_args.Get_String_Value("-o");
		driver.scale=parse_args.Get_Integer_Value("-s");
		driver.test=parse_args.Get_Integer_Value("-test");

		driver.Initialize();
		driver.Run();	
	}break;
	case 3:{//// Zalesak sphere test
		//// test 1 for particle levelset
		//// test 2 for standard levelset
		ParticleLevelsetAdvectionDriver<d> driver;
		driver.output_dir=parse_args.Get_String_Value("-o");
		driver.scale=100;
		driver.test=parse_args.Get_Integer_Value("-test");
		driver.last_frame=630; //// period=628
		driver.frame_rate=1;

		driver.Initialize();
		driver.Run();	
	}break;
	case 4:{//// fluid with particle levelset
		//// test 1 for droplet-on-tank scene
		//// test 2 for surface-tension scene
		FluidParticleLevelsetDriver<d> driver;
		driver.output_dir=parse_args.Get_String_Value("-o");
		driver.scale=parse_args.Get_Integer_Value("-s");
		driver.test=parse_args.Get_Integer_Value("-test");
		driver.fluid.use_particle=(parse_args.Get_Integer_Value("-use_pls")==1);

		driver.Initialize();
		driver.Run();	
	}break;
	}

}

#endif