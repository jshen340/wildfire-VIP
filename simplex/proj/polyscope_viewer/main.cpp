//////////////////////////////////////////////////////////////////////////
// Polyscope viewer
// Copyright (c) (2018-), Bo Zhu, Zhecheng Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "ParseArgs.h"
#include "AuxFunc.h"
#include "PolyscopeViewerDriver.h"
#include "Json.h"
#include "File.h"

int main(int argc,char* argv[])
{
    ////parse arguments
    ParseArgs parse_args;
	parse_args.Add_String_Argument("-o","./","data path");
	parse_args.Add_String_Argument("-oo","","rendering output");
	parse_args.Add_String_Argument("-m","base","viewer mode");
	parse_args.Add_String_Argument("-c","config.json","complex config settings");

    parse_args.Parse(argc,argv);

	std::string mode=parse_args.Get_String_Value("-m");

	std::shared_ptr<PolyscopeViewer> viewer=nullptr;

	////fluid viewers
	if(mode=="base"){viewer.reset(new PolyscopeViewer());}
	else if(mode=="fluid"){viewer.reset(new PolyscopeViewerEulerianFluid());}
	else{std::cout<<"Invalid viewer mode"<<std::endl;return 0;}

	//////////////////////////////////////////////////////////////////////////
	////These parameters need to be set before initialization
	viewer->output_dir=parse_args.Get_String_Value("-o");
	std::string render_dir=parse_args.Get_String_Value("-oo");
	if(render_dir=="")render_dir=parse_args.Get_String_Value("-o")+"/_images";
	viewer->render_dir=render_dir;
	std::string config=parse_args.Get_String_Value("-c");
	json config_json;
	if (File::File_Exists(config)) {
		std::ifstream in(config);
		config_json = json::parse(in);
	}

	//////////////////////////////////////////////////////////////////////////

	viewer->Initialize(config_json);

	viewer->Show();

	return 0;
}