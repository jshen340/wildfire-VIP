#include <iostream>
#include "ParseArgs.h"
#include "MeshToLevelSetDriver.h"
#ifndef __Main_cpp__
#define __Main_cpp__

int main(int argc,char* argv[])
{
	const int d = 3;

    ParseArgs parse_args;
    parse_args.Add_String_Argument("-o","output","output path");
	parse_args.Add_String_Argument("-input","input","intput mesh file path");
    parse_args.Add_Integer_Argument("-s",128,"resolution");
    parse_args.Add_Integer_Argument("-test",1,"test");
	parse_args.Add_Integer_Argument("-driver",1,"driver");

    parse_args.Parse(argc,argv);

    std::string output_dir=parse_args.Get_String_Value("-o");
	std::string input_path=parse_args.Get_String_Value("-input");
    const int scale=parse_args.Get_Integer_Value("-s");
	const int test=parse_args.Get_Integer_Value("-test");
	const int driver_id=parse_args.Get_Integer_Value("-driver");

	switch(driver_id){
		case 1: { // test MeshToLeveSet driver
			MeshToLevelSetDriver<d> driver;
			driver.scale = scale;
			driver.output_dir = output_dir;
			driver.input_path = input_path;
			driver.test = test;
			driver.Initialize();
			driver.Run();
		}break;
	}
}

#endif