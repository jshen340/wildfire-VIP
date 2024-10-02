//////////////////////////////////////////////////////////////////////////
// Diagonal Preconditioner for CPX
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "AuxFuncCPX.h"
#include "SparseMapping.h"

using namespace AuxFuncCPX;

template<class T>
class DiagonalPreconditioner: public LinearMapping {
public:
	T* diag_inv_dev = nullptr;
	int col_num = 0, row_num = 0;
	DiagonalPreconditioner(){}
	~DiagonalPreconditioner() {
		Global_Free(diag_inv_dev, DataHolder::DEVICE);
	}
	//TODO: support multi time Init() here
    void Init(SparseMapping<T>& A_dev);
    virtual int xDoF() { return col_num; }
    virtual int yDoF() { return row_num; }

    //input p, get Ap
	virtual void applyMapping(T* Ap_dev, T* p_dev);
};
