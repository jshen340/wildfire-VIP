//#####################################################################
// Main
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#include <iostream>
#include "ParseArgs.h"
#include "File.h"
#include "OptimizerDriver.h"

#ifndef __Main_cpp__
#define __Main_cpp__

int main(int argc,char* argv[])
{
    ////parse arguments
    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-o","output","output path");
    parse_args.Add_Integer_Argument("-s",8,"resolution");
	parse_args.Add_Integer_Argument("-test",1,"test");
    parse_args.Parse(argc,argv);

	std::string frame_dir=parse_args.Get_String_Value("-o");

	OptimizerTestDriver driver;
	driver.Initialize();
	driver.Run();
}

#endif