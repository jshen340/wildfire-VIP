//////////////////////////////////////////////////////////////////////////
// Sparse solver for CPX
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "SparseSolverCPX.h"

template<class T>
void SparseSolverCPX<T>::Init(const SparseMatrix<Scalar>& A_host, const int max_iter, const Scalar relative_tolerance)
{
	sparse_mapping.Update(A_host);
	pred.Init(sparse_mapping);
	cg = ConjugatedGradient();
	cg.preconditioner = &pred;
	cg.linear_mapping = &sparse_mapping;
	cg.Init(max_iter, relative_tolerance);
	checkCudaErrors(cudaGetLastError());
}

template<class T>
void SparseSolverCPX<T>::Solve(Scalar* x_host, Scalar const* b_host)
{
	cg.Solve(x_host, b_host);
}

template class SparseSolverCPX<Scalar>;