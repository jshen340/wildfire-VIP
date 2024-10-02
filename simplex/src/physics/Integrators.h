//////////////////////////////////////////////////////////////////////////
// Integrators for particle system
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "ParticleLambda.h"
#include "ImplicitShape.h"

//With V already known, integrate X.
//V is passed by reference because it may be changed in boundary treatment.
//We always want V to be consistent with the change of X.
template<int d>
class IntegratorX {
	Typedef_VectorD(d);
public:
	//IMPORTANT: make it parallel
	virtual void Integrate(Array<VectorD>& X, Array<VectorD>& V, const real dt) = 0;
};


//t_coeff=0: no-slip. t_coeff=1: free-slip
template<int d>
class IntegratorXImplicit : public IntegratorX<d> {
	Typedef_VectorD(d);
public:
	ImplicitShape<d> boundary;
	real t_coeff;
	real n_coeff;
	IntegratorXImplicit(const ImplicitShape<d> _boundary = ImplicitShape<d>(), const real _t_coeff = 1.0, const real _n_coeff = 0) {
		boundary = _boundary;
		t_coeff = _t_coeff;
		n_coeff = _n_coeff;
	}
	//see: A Variational Staggered Particle Framework for Incompressible Free-Surface Flows
	virtual void Integrate(Array<VectorD>& X, Array<VectorD>& V, const real dt);
};