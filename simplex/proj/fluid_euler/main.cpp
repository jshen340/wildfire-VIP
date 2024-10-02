//////////////////////////////////////////////////////////////////////////
// Main
// Copyright (c) (2018-), Bo Zhu, Yueyang Xianzang, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "ParseArgs.h"
#include "FluidEulerDriver.h"
#include "FluidEulerFreeSurfaceDriver.h"
#include "FluidPicFlipDriver.h"
#include "FluidViscousDriver.h"
#include "AdvectionDriver.h"
#include "FluidVorticityDriver.h"
#include "PoissonTestDriver.h"
#include "FluidAdvectionReflectionDriver.h"
#include "FluidEulerImmersedBoundaryDriver.h"
#include "FluidEulerTwoPhaseDriver.h"
//#include "CPXFunc.h"

#ifndef __Main_cpp__
#define __Main_cpp__

template<int d> void Run(ParseArgs& parse_args)
{

    std::string output_dir=parse_args.Get_String_Value("-o");
    const int scale=parse_args.Get_Integer_Value("-s");
	const int driver=parse_args.Get_Integer_Value("-driver");
	const int test=parse_args.Get_Integer_Value("-test");
	const int last_frame=parse_args.Get_Integer_Value("-lf");
	const real cfl=parse_args.Get_Double_Value("-cfl");
	const int frame_rate = parse_args.Get_Integer_Value("-fr");

	switch(driver){
	case 1:{
		try {
			FluidEulerDriver<d> driver;
			driver.scale = scale;
			driver.output_dir = output_dir;
			driver.test = test;
			driver.last_frame = last_frame;
			driver.frame_rate = frame_rate;
			driver.cfl = cfl;
			driver.Initialize();
			driver.Run();
		}
		catch (nlohmann::json::exception& e)
		{
			Info("json exception {}", e.what());
		}
	}break;
	case 2:{
		FluidEulerFreeSurfaceDriver<d> driver;
		driver.scale=scale;
		driver.output_dir=output_dir;
		driver.test=test;
		driver.last_frame=last_frame;
		driver.frame_rate = frame_rate;
		driver.cfl=cfl;
		driver.Initialize();
		driver.Run();		
	}break;
	case -2: {
#ifdef USE_CPX
		PoissonTestDriver Driver;
		Driver.grid_size = scale;
		Driver.test = test;
		Driver.Initialize();
		Driver.Run();
#endif
	}break;
	case 3:{
		FluidPicFlipDriver<d> driver;
		driver.scale=scale;
		driver.output_dir=output_dir;
		driver.test=test;
		driver.last_frame=last_frame;
		driver.frame_rate = frame_rate;
		driver.cfl=cfl;
		driver.Initialize();
		driver.Run();	
	}break;
	case 4:{
		FluidViscousDriver<d> driver;
		driver.scale=scale;
		driver.output_dir=output_dir;
		driver.test=test;
		driver.last_frame=last_frame;
		driver.frame_rate = frame_rate;
		driver.Initialize();
		driver.Run();	
	}break;	
	case 5:{
		AdvectionDriver<d> driver;
		driver.scale=scale;
		driver.output_dir=output_dir;
		driver.test=test;
		driver.last_frame=last_frame;
		driver.frame_rate = frame_rate;
		driver.Initialize();
		driver.Run();	
	}break;
	case 6:{
		////support 2D only
		FluidVorticityDriver<2> driver;
		driver.scale=scale;
		driver.output_dir=output_dir;
		driver.test=test;
		driver.last_frame=last_frame;
		driver.frame_rate = frame_rate;
		driver.Initialize();
		driver.Run();
	}break;
	case 7: {
		FluidAdvectionReflectionDriver<2> driver;
		driver.scale = scale;
		driver.output_dir = output_dir;
		driver.test = test;
		driver.last_frame = last_frame;
		driver.frame_rate = frame_rate;
		driver.cfl = cfl;
		driver.Initialize();
		driver.Run();
	}break;
	case 8: {
		FluidEulerImmersedBoundaryDriver<d> driver;
		driver.scale = scale;
		driver.output_dir = output_dir;
		driver.test = test;
		driver.last_frame = last_frame;
		driver.frame_rate = frame_rate;
		driver.cfl = cfl;
		driver.Initialize();
		driver.Run();
	}break;
	case 9: {
		FluidEulerTwoPhaseDriver<d> driver;
		driver.scale = scale;
		driver.output_dir = output_dir;
		driver.test = test;
		driver.last_frame = last_frame;
		driver.frame_rate = frame_rate;
		driver.cfl = cfl;
		driver.Initialize();
		driver.Run();
	}break;
	}
}

int main(int argc,char* argv[])
{
    ParseArgs parse_args;
    parse_args.Add_String_Argument("-o","output","output path");
	parse_args.Add_Integer_Argument("-s", 64, "resolution");
    parse_args.Add_Integer_Argument("-test",7,"test");
	parse_args.Add_Integer_Argument("-driver",2,"driver");
	parse_args.Add_Integer_Argument("-lf",100,"last frame");
	parse_args.Add_Integer_Argument("-fr",25,"frame rate");
	parse_args.Add_Double_Argument("-cfl",1.,"last frame");
	parse_args.Add_Integer_Argument("-d",2,"dimension");
    parse_args.Parse(argc,argv);

	const int d=parse_args.Get_Integer_Value("-d");
	
	if(d==2)Run<2>(parse_args);
	else if(d==3)Run<3>(parse_args);
}

#endif