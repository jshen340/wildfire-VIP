//////////////////////////////////////////////////////////////////////////
// Linear mapping represented by a sparse matrix
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "SparseMatrixCPX.h"
#include "AuxFuncCPX.h"
template<class T> class SparseMapping : public LinearMapping {
public:
	SparseMatrixCPX<Scalar> mat_dev;
	Scalar one = 1;
	Scalar zero = 0;
	cusparseHandle_t cusparseHandle = nullptr;

	SparseMapping() :mat_dev(DataHolder::DEVICE) {
		cusparseCreate(&cusparseHandle);
	}

	~SparseMapping() {
		if (cusparseHandle) cusparseDestroy(cusparseHandle);
	}

	void Update(const SparseMatrix<Scalar>& mat_host) {
		mat_dev.deepcopy(mat_host);
		cudaDeviceSynchronize();
		checkCudaErrors(cudaGetLastError());
	}

	//number of cols
	virtual int xDoF() {
		return mat_dev.cols();
	}

	//number of rows
	virtual int yDoF() {
		return mat_dev.rows();
	}

	//input p, get Ap
	virtual void applyMapping(Scalar* Ap_dev, Scalar* p_dev) {
		AuxFuncCPX::Global_Memset(Ap_dev, 0, yDoF(), DataHolder::DEVICE);
		AuxFuncCPX::Csrmv<Scalar>(cusparseHandle, &mat_dev, p_dev, Ap_dev, &one, &zero);
	}
};