//////////////////////////////////////////////////////////////////////////
// A modified version of fluid euler free surface to test the "lighthouse" scene
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "FluidEulerFreeSurface.h"

template<int d>
class FluidEulerLighthouse :public FluidEulerFreeSurface<d> {
public:
	Typedef_VectorDii(d);
	using Base = FluidEulerFreeSurface<d>;
	real 
};