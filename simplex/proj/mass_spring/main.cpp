//#####################################################################
// Main
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com, Yankai Mao
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#include <iostream>
#include "ParseArgs.h"
#include "MassSpringDriver.h"
#include "ClothDriver.h"

#ifndef __Main_cpp__
#define __Main_cpp__

//#define ROOT_PATH2 123
int main(int argc,char* argv[])
{
    const int d=3;

    ParseArgs parse_args;
    parse_args.Add_String_Argument("-o","output","output path");
    parse_args.Add_Integer_Argument("-s",1,"resolution");
    parse_args.Add_Integer_Argument("-test",1,"test");
	parse_args.Add_Integer_Argument("-driver",1,"driver");
	parse_args.Add_Integer_Argument("-lf",200,"last frame");
    parse_args.Parse(argc,argv);

    std::string output_dir=parse_args.Get_String_Value("-o");
    const int scale=parse_args.Get_Integer_Value("-s");
	const int driver=parse_args.Get_Integer_Value("-driver");
	const int test=parse_args.Get_Integer_Value("-test");
	const int last_frame=parse_args.Get_Integer_Value("-lf");

	switch(driver){
	case 1:{
		MassSpringDriver<d> driver;
		driver.scale=scale;
		driver.output_dir=output_dir;
		driver.test=test;
		driver.last_frame=last_frame;
		driver.soft_body.use_newmark=true;
		driver.Initialize();
		driver.Run();	
	}break;
	case 2:{
		ClothDriver driver;
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