#include <iostream>
#include "ParseArgs.h"
#include "FluidMicroDriver.h"

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

	switch(driver){
	case 1:{
		FluidMicroDriver<d> driver;
		driver.scale=scale;
		driver.output_dir=output_dir;
		driver.test=test;
		driver.last_frame=last_frame;
		driver.cfl=cfl;
		driver.Initialize();
		driver.Run();	
	}break;
	}
}

int main(int argc,char* argv[])
{
    ParseArgs parse_args;
    parse_args.Add_String_Argument("-o","output","output path");
    parse_args.Add_Integer_Argument("-s",64,"resolution");
    parse_args.Add_Integer_Argument("-test",-2,"test");
	parse_args.Add_Integer_Argument("-driver",1,"driver");
	parse_args.Add_Integer_Argument("-lf",200,"last frame");
	parse_args.Add_Double_Argument("-cfl",1.,"last frame");
	parse_args.Add_Integer_Argument("-d",2,"dimension");
    parse_args.Parse(argc,argv);

	const int d=parse_args.Get_Integer_Value("-d");
	
	if(d==2)Run<2>(parse_args);
	else if(d==3)Run<3>(parse_args);
}

#endif