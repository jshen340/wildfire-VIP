//#####################################################################
// Poisson main
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#include <iostream>
#include "ParseArgs.h"
#include "Common.h"
#include "PoissonDriver.h"
#include "PoissonIrregularDriver.h"

#ifndef __Main_cpp__
#define __Main_cpp__

int main(int argc,char* argv[])
{
	const int d=2;
    ParseArgs parse_args;
    parse_args.Add_String_Argument("-o","output","output path");
    parse_args.Add_Integer_Argument("-s",32,"resolution");
	parse_args.Add_Integer_Argument("-test",2,"test source");
	parse_args.Add_Integer_Argument("-driver",1,"driver id");
	parse_args.Add_Integer_Argument("-boundary",1,"irregular boundary type");
	parse_args.Add_Integer_Argument("-boundary_type",1,"D or N bc type");
    parse_args.Parse(argc,argv);
	
	int driver=parse_args.Get_Integer_Value("-driver");

	switch(driver){
	case 1:{
		PoissonDriver<d> driver;
		driver.output_directory=parse_args.Get_String_Value("-o");
		driver.scale=parse_args.Get_Integer_Value("-s");
		driver.test=parse_args.Get_Integer_Value("-test");
		driver.Initialize();
		driver.Run();	
	}break;
	case 2:{
		PoissonIrregularDriver<d> driver;
		driver.output_directory=parse_args.Get_String_Value("-o");
		driver.scale=parse_args.Get_Integer_Value("-s");
		driver.test=parse_args.Get_Integer_Value("-test");
		driver.boundary=parse_args.Get_Integer_Value("-boundary"); 
		driver.boundary_type=parse_args.Get_Integer_Value("-boundary_type");
		driver.Initialize();
		driver.Run();
	}break;
	}

}

#endif