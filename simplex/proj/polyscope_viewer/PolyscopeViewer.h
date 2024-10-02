//////////////////////////////////////////////////////////////////////////
// Polyscope viewer
// Copyright (c) (2018-), Bo Zhu, Zhecheng Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __PolyscopeViewer_h__
#define __PolyscopeViewer_h__
#include "Common.h"
#include "Json.h"
#include "PolyscopeGrid.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

class PolyscopeViewer
{
public:
	json j_fields;
	std::string output_dir,render_dir;
	int first_frame=0,last_frame=-1,frame=0,loaded_frame=-1;
	bool play=false;
	bool render=false;
	bool verbose=true;

	//////////////////////////////////////////////////////////////////////////
	////Initialization and run
	virtual void Initialize(const json& j);
	virtual void Register_Common_Data();
	virtual void Register_Data();
	virtual void Update_Data();
	virtual void Register_UI();
	virtual void Run();
	virtual void Render();
	virtual void Show();

	//////////////////////////////////////////////////////////////////////////
	////Register Polyscope structures
	virtual polyscope::PointCloud* Register_Point_Cloud(std::string name, bool enabled, glm::vec3 color);
	virtual polyscope::PointCloud* Update_Point_Cloud(std::string name);
	virtual PolyscopeGrid* Register_Grid(std::string name, bool enabled);
	virtual void Register_Grid_Cell_Scalar_Field(std::string grid_name, std::string name, bool enabled, std::string color_map);
	virtual void Update_Grid_Cell_Scalar_Field(std::string grid_name, std::string name);
	virtual void Register_Grid_Cell_Vector_Field(std::string grid_name, std::string name, bool enabled, glm::vec3 color);
	virtual void Update_Grid_Cell_Vector_Field(std::string grid_name, std::string name);
	virtual polyscope::CurveNetwork* Register_Surface_Mesh_2D(std::string name, bool enabled, glm::vec3 color);
	virtual polyscope::CurveNetwork* Update_Surface_Mesh_2D(std::string name);
	virtual polyscope::SurfaceMesh* Register_Surface_Mesh(std::string name, bool enabled, glm::vec3 color);
	virtual polyscope::SurfaceMesh* Update_Surface_Mesh(std::string name);

	//////////////////////////////////////////////////////////////////////////
	////Set object properties
	template<class T_OBJECT> void Set_Visibility(T_OBJECT* obj, const bool visible = true) {
		if (obj == nullptr)return;
		obj->setEnabled(visible);
	}

	template<class T_OBJECT> void Set_Transparancy(T_OBJECT* obj, const float alpha = 1.f) {
		if (obj == nullptr)return;
		obj->setTransparancy(alpha);
	}

	//////////////////////////////////////////////////////////////////////////
	////Helper functions
	void Extract_Points(Array<float>& gl_array, Array<Vector3f>& points);
	void Update_Last_Frame();
};

#endif
