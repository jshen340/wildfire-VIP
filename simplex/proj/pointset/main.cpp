#include <iostream>
#include "ParseArgs.h"
#include "PointSetDriver.h"
#include "LeastSquaresDriver.h"
#include "MlsOptimizerDriver.h"

#ifndef __Main_cpp__
#define __Main_cpp__

int main(int argc,char* argv[])
{
    const int d=3;	//

    ParseArgs parse_args;
    parse_args.Add_String_Argument("-o","output","output path");
    parse_args.Add_Integer_Argument("-s",64,"resolution");
    parse_args.Add_Integer_Argument("-test",6,"test");
	parse_args.Add_Integer_Argument("-driver",2,"driver");
	parse_args.Add_Integer_Argument("-lf",200,"last frame");
	parse_args.Add_Double_Argument("-cfl",1,"last frame");
    parse_args.Parse(argc,argv);

    std::string output_dir=parse_args.Get_String_Value("-o");
    const int scale=parse_args.Get_Integer_Value("-s");
	const int driver=parse_args.Get_Integer_Value("-driver");
	const int test=parse_args.Get_Integer_Value("-test");
	const int last_frame=parse_args.Get_Integer_Value("-lf");
	const real cfl=parse_args.Get_Double_Value("-cfl");

	switch(driver){
	case 1:{
		PointSetDriver<d> driver;
		driver.scale=scale;
		driver.output_dir=output_dir;
		driver.test=test;
		driver.last_frame=last_frame;
		driver.cfl=cfl;
		driver.Initialize();
		driver.Run();	
	}break;
	case 2:{
		LeastSquaresDriver<d> driver;
		driver.scale=scale;
		driver.output_dir=output_dir;
		driver.test=test;
		driver.last_frame=last_frame;
		driver.cfl=cfl;
		driver.Initialize();
		driver.Run();
	}break;
	case 3:{
		MlsOptimizerDriver<d> driver;
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

#endif