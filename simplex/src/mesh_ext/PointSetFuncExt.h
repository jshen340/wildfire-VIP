#pragma once

#include "Mesh.h"
#include "MeshFunc.h"
#include "AuxFunc.h"
#include "Constants.h"
#include "GeometryParticles.h"
#include "RandomNumber.h"
#include "MacGrid.h"
#include "GeometryPrimitives.h"
#include "RenderFunc.h"

#include "MeshFuncExt.h"

namespace PointSetFunc
{
	using namespace AuxFunc;
    //////////////////////////////////////////////////////////////////////////
	////3D points initialization
	real Initialize_Sphere_Points(const Vector3& c, const real R, const int sub, GeometryParticles<3>& particles, const real randomness = 0.);
    real Initialize_Half_Sphere_Points(const Vector3& c, const real R, const int sub, GeometryParticles<3>& particles, std::vector<int>* is_boundary);
    real Add_Sphere_Points(const Vector3& c, const real R, const int sub, GeometryParticles<3>& particles);
    real Initialize_Circle_Points(const Vector3& c, const real R, const Vector3& normal, const real dx, GeometryParticles<3>& particles);
    real Initialize_Ellipsoid_Points(const Vector3& C, const real R, const int sub, const real a, const real b, const real c, GeometryParticles<3>& particles);
}