//////////////////////////////////////////////////////////////////////////
// Diagonal Preconditioner for CPX
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "DiagonalPreconditioner.h"
#include "gpuUtils.h"

template<class T>
void DiagonalPreconditioner<T>::Init(SparseMapping<T>& A_dev) {
    col_num = A_dev.xDoF();
    row_num = A_dev.yDoF();
    Assert(col_num == row_num, "DiagonalPreconditioner::Init column number doesn't equal to row number");
    Assert(col_num > 0, "DiagonalPreconditioner::Init can't solve an empty system");
    diag_inv_dev = Global_Malloc<T>(col_num, DataHolder::DEVICE);
    
    const int* ptr_dev = A_dev.mat_dev.outIndexPtr();
    const int* col_dev = A_dev.mat_dev.innerIndexPtr();
    const T* val_dev = A_dev.mat_dev.valuePtr();
    //note: if you write something like A_dev.mat_dev.ptr in the lambda function, it will fail
    //seems __device__ lambda can't properly capture object member
    auto f = [ptr_dev,col_dev,val_dev]__device__(T & inv_i, const int i) {
        //row of i
        int row_start = ptr_dev[i], row_end = ptr_dev[i + 1];
        bool found = false;
        for (int it = row_start;it < row_end;it++) {
            int col_num = col_dev[it];
            if (col_num == i) {
                T val = val_dev[it];
                if (val != 0) {
                    found = true;
                    inv_i = 1.0 / val;
                }
                break;
            }
            else if (col_num > i) break;
        }
        if (!found) {
            inv_i = 1.0;
        }
    };
    cwise_mapping_with_idx_wrapper(diag_inv_dev, f, row_num);
    checkCudaErrors(cudaGetLastError());
}

template<class T>
void DiagonalPreconditioner<T>::applyMapping(T* Ap_dev, T* p_dev) {
    Scalar* inv_dev = this->diag_inv_dev;
    auto f = [inv_dev, p_dev]__device__(T & Ap_i, const int i) {
        Ap_i = p_dev[i] * inv_dev[i];
    };
    cwise_mapping_with_idx_wrapper(Ap_dev, f, row_num);
}

template class DiagonalPreconditioner<Scalar>;
