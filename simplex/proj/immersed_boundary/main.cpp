//////////////////////////////////////////////////////////////////////////
// Mesh To Level Set Driver
// Copyright (c) (2021-), Fan Feng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "ParseArgs.h"
#include "ImmersedBoundaryDriver.h"
#ifndef __Main_cpp__
#define __Main_cpp__

template<int d> void Run_Program(ParseArgs& parse_args)
{
	std::string output_dir = parse_args.Get_String_Value("-o");
	std::string input_path = parse_args.Get_String_Value("-input");
	const int scale = parse_args.Get_Integer_Value("-s");
	const int test = parse_args.Get_Integer_Value("-test");
	const int driver_id = parse_args.Get_Integer_Value("-driver");
	const int last_frame = parse_args.Get_Integer_Value("-lf");
	const real cfl = parse_args.Get_Double_Value("-cfl");

	switch (driver_id) {
	case 1: {
		ImmersedBoundaryDriver<d> driver;
		driver.scale = scale;
		driver.output_dir = output_dir;
		driver.test = test;
		driver.last_frame = last_frame;
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
	parse_args.Add_String_Argument("-input","input","intput mesh file path");
    parse_args.Add_Integer_Argument("-s",128,"resolution");
    parse_args.Add_Integer_Argument("-test",1,"test");
	parse_args.Add_Integer_Argument("-driver",1,"driver");
	parse_args.Add_Integer_Argument("-d", 2, "dimension");
	parse_args.Add_Integer_Argument("-lf", 100, "last frame");
	parse_args.Add_Double_Argument("-cfl", 1, "cfl");

    parse_args.Parse(argc,argv);

	std::string output_dir = parse_args.Get_String_Value("-o");
	if (!File::Directory_Exists(output_dir.c_str()))File::Create_Directory(output_dir);
	std::string command_file = output_dir + "/command.txt";
	File::Write_Text_To_File(command_file, *argv);
	while (--argc) { File::Append_Text_To_File(command_file, std::string(" ") + *(++argv)); }

	int d = parse_args.Get_Integer_Value("-d");
	if (d == 2)Run_Program<2>(parse_args);
	else Run_Program<3>(parse_args);
}

#endif