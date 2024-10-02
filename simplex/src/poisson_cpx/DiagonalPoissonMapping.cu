//////////////////////////////////////////////////////////////////////////
// Linear mapping with cell itself and neighbors on grid
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "DiagonalPoissonMapping.h"
#include "gpuUtils.h"
#include "form01mapping.h"
#include "AuxFuncCPX.h"
#include <fmt/os.h>
#include <cuda.h>


template<class T, int d>
void DiagonalPoissonMapping<T, d>::Init(const VectorDi _N)
{
	N = _N;
	dof = N.prod();

	descr.init(N);
	poisson_mapping.init(&descr);

	additional_diag_dev = AuxFuncCPX::Global_Malloc<T>(dof, DataHolder::DEVICE);
	additional_diag_host = AuxFuncCPX::Global_Malloc<T>(dof, DataHolder::HOST);
}

template<class T, int d>
void DiagonalPoissonMapping<T, d>::applyMapping(Scalar* Ap_dev, Scalar* p_dev) {
	poisson_mapping.applyMapping(Ap_dev, p_dev);
	auto fix_diag_func = [=]__device__(Scalar & alpha_i, bool fixed_i) { if (fixed_i) alpha_i = 0; };
	auto add_diagonal_func = [=] __device__(Scalar & Ap_i, Scalar p_i, Scalar alpha_i) { Ap_i += p_i * alpha_i; };
	cwise_mapping_wrapper(additional_diag_dev, descr.d_fixed, fix_diag_func, dof);
	cwise_mapping_wrapper(Ap_dev, p_dev, additional_diag_dev, add_diagonal_func, dof);
}

template<class T, int d>
void DiagonalPoissonMapping<T, d>::To_Device(void)
{
	AuxFuncCPX::Global_Copy_Array(additional_diag_dev, additional_diag_host, dof, DataHolder::DEVICE, DataHolder::HOST);
	descr.toDevice();
	descr.finish();
}


template class DiagonalPoissonMapping<Scalar, 2>;
template class DiagonalPoissonMapping<Scalar, 3>;

// for blockDim = (8, 8)
// iterate through cell
// set cell(i,j) to the sum of all its neighboring faces
// seems that it behaves like cell(i,j)=0 for outside the boundary
__global__ void Sum_Faces_To_Cell_Kernel2(grid2D grid, Scalar* face_x, Scalar* face_y, Scalar* cell_output)
{
	//const int nbx = gridDim.x;
	//const int nby = gridDim.y;

	//index of block in all 8*8 blocks
	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	//index of thread in the particular 8*8 block
	const int idx = threadIdx.x;
	const int idy = threadIdx.y;

	//8*9 or 9*8 face data
	__shared__ Scalar shared_face_data[72];

	Scalar sum = 0;
	{
		// 9x8 face_x data load
		// shared_face_data(idx,idy)=
		shared_face_data[idy * 8 + idx] = face_x[grid.face_ind(bx * 8 + idy, by * 8 + idx, 0)];
		if (idy == 0) shared_face_data[64 + idx] = face_x[grid.face_ind((bx + 1) * 8 + idy, by * 8 + idx, 0)];
		__syncthreads();

		// left x-axis faces
		sum += shared_face_data[idx * 8 + idy];

		// right x-axis faces
		sum += shared_face_data[(idx + 1) * 8 + idy];
	}
	__syncthreads();

	{
		// 8x9 face_y data load
		shared_face_data[idy * 8 + idx] = face_y[grid.face_ind(bx * 8 + idx, by * 8 + idy, 1)];
		if (idy == 0) shared_face_data[64 + idx] = face_y[grid.face_ind(bx * 8 + idx, (by + 1) * 8 + idy, 1)];
		__syncthreads();

		// down y-axis faces
		sum += shared_face_data[idy * 8 + idx];

		// up y-axis faces
		sum += shared_face_data[(idy + 1) * 8 + idx];
	}

	cell_output[grid.cell_ind(bx * 8 + idx, by * 8 + idy)] = sum;
}

// for blockDim = (4, 4, 4)
// iterate through cell
// set cell(i, j, k) to the sum of all its neighboring faces
__global__ void Sum_Faces_To_Cell_Kernel3(grid3D grid, Scalar* face_x, Scalar* face_y, Scalar* face_z, Scalar* cell_output)
{
	//const int nbx = gridDim.x;
	//const int nby = gridDim.y;
	//const int nbz = gridDim.z;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;
	const int bz = blockIdx.z;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;
	const int idz = threadIdx.z;

	//5*4*4 or 4*5*4 or 4*4*5 face data
	__shared__ Scalar shared_face_data[80];

	Scalar sum = 0;
	{
		// 5x4x4 face_x data load
		//shared_face_data(tx,ty,tz)=face_x_block(tz,ty,tx)
		shared_face_data[(idz * 4 + idy) * 4 + idx] = face_x[grid.face_ind(bx * 4 + idz, by * 4 + idx, bz * 4 + idy, 0)];
		if (idz == 0) shared_face_data[(4 * 4 + idy) * 4 + idx] = face_x[grid.face_ind(bx * 4 + 4, by * 4 + idx, bz * 4 + idy, 0)];
		__syncthreads();

		// left x-axis faces
		sum += shared_face_data[(idx * 4 + idz) * 4 + idy];

		// right x-axis faces
		sum += shared_face_data[((idx + 1) * 4 + idz) * 4 + idy];
	}
	__syncthreads();

	{
		// 4x5x4 face_y data load
		shared_face_data[(idz * 4 + idy) * 4 + idx] = face_y[grid.face_ind(bx * 4 + idx, by * 4 + idz, bz * 4 + idy, 1)];
		if (idz == 0) shared_face_data[(4 * 4 + idy) * 4 + idx] = face_y[grid.face_ind(bx * 4 + idx, by * 4 + 4, bz * 4 + idy, 1)];
		__syncthreads();

		// down y-axis faces
		sum += shared_face_data[(idy * 4 + idz) * 4 + idx];

		// up y-axis faces
		sum += shared_face_data[((idy + 1) * 4 + idz) * 4 + idx];
	}
	__syncthreads();

	{
		// 4x4x5 face_y data load
		shared_face_data[(idz * 4 + idy) * 4 + idx] = face_z[grid.face_ind(bx * 4 + idx, by * 4 + idy, bz * 4 + idz, 2)];
		if (idz == 0) shared_face_data[(4 * 4 + idy) * 4 + idx] = face_z[grid.face_ind(bx * 4 + idx, by * 4 + idy, bz * 4 + 4, 2)];
		__syncthreads();

		// front z-axis faces
		sum += shared_face_data[(idz * 4 + idy) * 4 + idx];

		// back y-axis faces
		sum += shared_face_data[((idz + 1) * 4 + idy) * 4 + idx];
	}
	__syncthreads();

	cell_output[grid.cell_ind(bx * 4 + idx, by * 4 + idy, bz * 4 + idz)] = sum;
}

template<class T, int d>
void DiagonalPoissonPreconditioner<T, d>::Init(const int _dof)
{
	dof = _dof;
	AuxFuncCPX::Global_Free(diag_inv_dev, DataHolder::DEVICE);
	diag_inv_dev = AuxFuncCPX::Global_Malloc<T>(dof, DataHolder::DEVICE);
}

template<class T, int d>
void DiagonalPoissonPreconditioner<T, d>::Compute(DiagonalPoissonMapping<T, d>& A_dev)
{
	if (A_dev.xDoF() != dof) Init(A_dev.xDoF());
	if constexpr (d == 2) {
		grid2D grid = A_dev.descr.grid;
		Sum_Faces_To_Cell_Kernel2 << < dim3((grid.Nx >> 3), (grid.Ny >> 3)), dim3(8, 8) >> > (
			grid,
			A_dev.descr.d_vol,
			A_dev.descr.d_vol + grid.face_size(0),
			diag_inv_dev
			);
	}
	else if constexpr (d == 3) {
		grid3D grid = A_dev.descr.grid;
		//grid3DOperator::D2Mapping(grid, A_dev.descr.d_vol,
		//	A_dev.descr.d_vol + grid.face_size(0),
		//	A_dev.descr.d_vol + grid.face_size(1),
		//	diag_inv_dev);
		Sum_Faces_To_Cell_Kernel3 << <dim3((grid.Nx >> 2), (grid.Ny >> 2), (grid.Nz >> 2)), dim3(4, 4, 4) >> > (
			grid,
			A_dev.descr.d_vol,
			A_dev.descr.d_vol + grid.face_size(0),
			A_dev.descr.d_vol + grid.face_size(0) + grid.face_size(1),
			diag_inv_dev
			);
	}
	auto plus_diag_extra = [=]__device__(T & inv_dev_i, T diag_extra_i) {
		inv_dev_i += diag_extra_i;
	};
	cwise_mapping_wrapper(diag_inv_dev, A_dev.additional_diag_dev, plus_diag_extra, dof);
	auto calc_inv_with_fixed = [=]__device__(T & inv_dev_i, bool fixed_dev_i) {
		if (fixed_dev_i || inv_dev_i == 0) {
			inv_dev_i = 1.0;
		}
		else inv_dev_i = 1.0 / inv_dev_i;
	};
	cwise_mapping_wrapper(diag_inv_dev, A_dev.descr.d_fixed, calc_inv_with_fixed, dof);
}

template<class T, int d>
void DiagonalPoissonPreconditioner<T, d>::applyMapping(T* Ap_dev, T* p_dev)
{
	Scalar* inv_dev = this->diag_inv_dev;
	auto f = [inv_dev, p_dev]__device__(T & Ap_i, const int i) {
		Ap_i = p_dev[i] * inv_dev[i];
	};
	cwise_mapping_with_idx_wrapper(Ap_dev, f, dof);
}

template class DiagonalPoissonPreconditioner<Scalar, 2>;
template class DiagonalPoissonPreconditioner<Scalar, 3>;