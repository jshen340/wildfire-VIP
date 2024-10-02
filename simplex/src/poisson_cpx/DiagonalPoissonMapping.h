//////////////////////////////////////////////////////////////////////////
// Linear mapping with cell itself and neighbors on grid
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "form01mapping.h"
#include "grid3D.h"
#include "PoissonDescriptor.h"
#include "TypeFunc.h"
#include "PoissonMapping.h"
#include "PoissonMapping3D.h"
#include "MacGrid.h"
template<class T, int d>
class DiagonalPoissonMapping : public LinearMapping {
	Typedef_VectorDii(d);
public:
	VectorDi N;
	int dof;//number of cells

	PoissonDescriptor<d> descr;
	using FixedMapping = typename If<d == 2, PoissonMappingFixed, PoissonMapping3DFixed >::Type;
	FixedMapping poisson_mapping;
	
	T* additional_diag_dev = nullptr, * additional_diag_host = nullptr;

	~DiagonalPoissonMapping() {
		AuxFuncCPX::Global_Free(additional_diag_dev, DataHolder::DEVICE);
		AuxFuncCPX::Global_Free(additional_diag_host, DataHolder::HOST);
	}

	void Init(const VectorDi _N);

	virtual int xDoF() { return dof; }

	virtual int yDoF() { return dof; }

	//input p, get Ap
	virtual void applyMapping(Scalar* Ap, Scalar* p);

	void To_Device(void);
};

//it's a diagonal preconditioner
template<class T, int d>
class DiagonalPoissonPreconditioner :public LinearMapping {
public:
	T* diag_inv_dev = nullptr;
	int dof = 0;
	DiagonalPoissonPreconditioner() {}
	~DiagonalPoissonPreconditioner() {
		AuxFuncCPX::Global_Free(diag_inv_dev, DataHolder::DEVICE);
	}
	void Init(const int _dof);
	void Compute(DiagonalPoissonMapping<T, d>& A_dev);
	virtual int xDoF() { return dof; }
	virtual int yDoF() { return dof; }

	virtual void applyMapping(T* Ap_dev, T* p_dev);
};