//////////////////////////////////////////////////////////////////////////
// Sparse solver for CPX
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "SparseMapping.h"
#include "ConjugatedGradient.h"
#include "DiagonalPreconditioner.h"

template<class T>
class SparseSolverCPX {
public:
	SparseMapping<T> sparse_mapping;
	DiagonalPreconditioner<T> pred;
	ConjugatedGradient cg;

	void Init(const SparseMatrix<Scalar>& A_host, const int max_iter = 50, const Scalar relative_tolerance = 1e-5);
	void Solve(Scalar* x_host, Scalar const* b_host);
};

