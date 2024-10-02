//#####################################################################
// Main
// Copyright (c) (2018-), Xiangxin Kong, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//#####################################################################
#include <iostream>
#include "ParseArgs.h"
#include "FluidSPHDriver.h"
#include "Params.h"
#include "SPHDriver.h"

#ifndef __Main_cpp__
#define __Main_cpp__

void Set_Threads(ParseArgs& parse_args) {
	int number_threads = parse_args.Get_Integer_Value("-tnum");
	omp_set_num_threads(number_threads);
	int max_threads = omp_get_max_threads();
	std::cout << "#     Set " << number_threads << " threads, run with " << max_threads << " cores\n";
}

int main(int argc,char* argv[])
{
	const int d = 2;

    ParseArgs parse_args;
    parse_args.Add_String_Argument("-o","output","output path");
    parse_args.Add_Integer_Argument("-s",1000,"resolution");
    parse_args.Add_Integer_Argument("-test",3,"test");
	parse_args.Add_Integer_Argument("-driver",1,"driver");
	parse_args.Add_Integer_Argument("-lf",1,"last frame");
	parse_args.Add_Integer_Argument("-algo", 1, "algorithm type");
	parse_args.Add_Integer_Argument("-tnum", 1, "number of threads");
    parse_args.Parse(argc,argv);
	
	Set_Threads(parse_args);

    std::string output_dir=parse_args.Get_String_Value("-o");
    const int scale=parse_args.Get_Integer_Value("-s");
	const int driver=parse_args.Get_Integer_Value("-driver");
	const int test=parse_args.Get_Integer_Value("-test");
	const int last_frame=parse_args.Get_Integer_Value("-lf");
	switch(driver){
	case 1: {
		SPHDriver<d> driver;
		driver.scale = scale;
		driver.output_dir = output_dir;
		driver.test = test;
		driver.last_frame = last_frame;
		driver.algo_type = parse_args.Get_Integer_Value("-algo");
		driver.Initialize();
		driver.Run();
	}break;
	case 2: {
		FluidSPHDriver<d> driver;
		driver.scale=scale;
		driver.output_dir=output_dir;
		driver.test=test;
		driver.last_frame=last_frame;
		driver.Initialize();
		driver.Run();
	}break;
	}
}

#endif