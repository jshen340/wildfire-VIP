//////////////////////////////////////////////////////////////////////////
// Auxiliary Function CUDA
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "SparseMatrixCPX.h"
//#include "gpuUtils.h"
#include "TypeFunc.h"
#include "form01mapping.h"
#include "grid3D.h"
#include "MacGrid.h"
#include <cublas_v2.h>

namespace AuxFuncCPX
{
	//////////////////////////////////////////////////////////////////////////
	////Basic APIs
	
	////CPX use 8*8 block in 2D and 4*4*4 block in 3D
	int CUDA_Block_Size(const int d);

	////Round up vector to a multiple of bn
	template<int d> Vector<int, d> Round_Up_To_Align(Vector<int, d> v, int bn) {
		for (int i = 0;i < d;i++) {
			//round to the nearest multiple of block_size
			v[i] = ((v[i] + bn - 1) / bn) * bn;
		}
		return v;
	}
	 
	////Data holder calculation
	DataHolder Copy_Data_Holder(const DataHolder& to, const DataHolder& from);

	////CUDA type selector
	template<class T> cudaDataType_t Cuda_Real_Type(void);

	////Memory allocation, communication between CPU and GPU
	//"Global" comes from __global__ in CUDA, for it can handle both device and host situations
	template<class T> T* Global_Malloc(int n, const DataHolder &side);
	template<class T> T* Global_Free(T*& ptr, const DataHolder &side);//Returns nullptr. When receiving nullptr, it does nothing
	template<class T> void Global_Copy_Array(T*& ptr_to, T*& ptr_from, int n, const DataHolder& to_side, const DataHolder& from_side);
	template<class T> T* Global_Malloc_And_Copy_Array(T*& ptr_from, int n, const DataHolder& to_side, const DataHolder& from_side);
	template<class T> void Global_Realloc_And_Copy_Array(T*& ptr_to, T*& ptr_from, int n, const DataHolder& to_side, const DataHolder& from_side);
	template<class T> void Global_Realloc_And_Copy(T*& ptr_to, T*& ptr_from, const DataHolder& to_side, const DataHolder& from_side);//only 1 element
	template<class T> void Global_Memset(T* ptr, int v, int n, const DataHolder& side);
	//////////////////////////////////////////////////////////////////////////
	////Linear algebra operations on device
	////sparse matrix-vector multiplication
	template<class T> void Csrmv(cusparseHandle_t handle, SparseMatrixCPX<T>* A_dev, T* x_dev, T* y_dev, T* alpha, T* beta);//T={double, float}, y=alpha*A*x+beta*y
	//////matrix-vector multiplication
	//void Mv(SparseMatrixCPX<double>* A_dev,double* x_dev,double* y_dev);
	//void Mv(SparseMatrixCPX<float>* A_dev,float* x_dev,float* y_dev);
	//////matrix-matrix multiplication
	//template<class T> void SpGEMM(SparseMatrixCPX<T>* A_dev, SparseMatrixCPX<T>* B_dev, SparseMatrixCPX<T>* C_dev);//T={double, float}, C=A*B
	//void Mm(SparseMatrixCPX<double>* A_dev, SparseMatrixCPX<double>* B_dev, SparseMatrixCPX<double>* C_dev);
	//void Mm(SparseMatrixCPX<float>* A_dev, SparseMatrixCPX<float>* B_dev, SparseMatrixCPX<float>* C_dev);
	//////copy
	//void Copy(double* to,double* from,int n);
	//void Copy(float* to,float* from,int n);

	template<class T> cudaDataType_t Cuda_Real_Type(void) 
	{
		int siz = sizeof(T);
		if (siz == 4) { return CUDA_R_32F; }
		else if (siz == 8) { return CUDA_R_64F; }
		else { std::cerr << "[Error] AuxFuncCuda::Cuda_Type: Unknown data type\n"; return cudaDataType_t(); }
	}

	// In the domain of fill_grid, fill CPX grid saved in linear form F_cpx_host with function F
	// Note that fill_grid can be smaller than grid_cpx. This is deliberately designed, for CPX will pad length to multiple of 4 or 8
	// Fixed the compiling issue for CUDA 11.1: changing std::function<T(const Vector<int, d>)> to F
	template<class T, int d, class F> 
	void Fill_CPX_Grid(const Grid<d> &fill_grid, F f,typename If<d == 2, grid2D, grid3D >::Type grid_cpx, T* F_cpx_host) {
		int cell_num = fill_grid.Number_Of_Cells();
#pragma omp parallel for
		for (int i = 0; i < cell_num; i++) {
			Vector<int, d> cell = fill_grid.Cell_Coord(i);
			int cell_ind = -1;
			if constexpr (d == 2) cell_ind = grid_cpx.cell_ind(cell[0], cell[1]);
			else if constexpr (d == 3) cell_ind = grid_cpx.cell_ind(cell[0], cell[1], cell[2]);
			F_cpx_host[cell_ind] = f(cell);
		}
	}

	template<class T, int d>
	void Fill_CPX_Face(const MacGrid<d>& fill_grid, std::function<T(const int i, const Vector<int, d>&)> F,typename If<d == 2, grid2D, grid3D>::Type grid_cpx, T* F_cpx_host) {
		int offset = 0;
		for (int axis = 0; axis < d; axis++) {
			int face_num = fill_grid.Number_Of_Faces(axis);
#pragma omp parallel for
			for (int i = 0; i < face_num; i++) {
				Vector<int, d> face = fill_grid.Face_Coord(axis, i);
				int face_ind = -1;
				if constexpr (d == 2) { face_ind = offset + grid_cpx.face_ind(face[0], face[1], axis); }
				else if constexpr (d == 3) { face_ind = offset + grid_cpx.face_ind(face[0], face[1], face[2], axis); }
				F_cpx_host[face_ind] = F(axis, face);
			}
			offset += grid_cpx.face_size(axis);
		}
	}
};
