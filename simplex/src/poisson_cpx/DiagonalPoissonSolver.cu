#include "DiagonalPoissonSolver.h"
#include "SimplexSupport.h"
#include "form01mapping.h"
#include "gpuUtils.h"
#include "AuxFunc.h"

template<class T, int d>
void DiagonalPoissonSolver<T, d>::Init(const MacGrid<d>& mac_grid, const int max_iter, const Scalar relative_tolerance)
{
	physical_grid = mac_grid;
	//grid_size is aligned, so may be larger than physical_grid
	VectorDi grid_size =
		AuxFuncCPX::Round_Up_To_Align<d>(
			physical_grid.grid.cell_counts,
			AuxFuncCPX::CUDA_Block_Size(d)
			);
	x_field_host.Resize(grid_size);

	//Initialize linear mapping
	diagonal_mapping.Init(grid_size);
	//set the data beyond padding area
#pragma omp parallel for
	for (int i = 0;i < diagonal_mapping.descr.fsize;i++) diagonal_mapping.descr.h_vol[i] = (Scalar)0;
#pragma omp parallel for
	for (int i = 0;i < diagonal_mapping.descr.size;i++) diagonal_mapping.descr.h_fixed[i] = true;
	x_linear_host = AuxFuncCPX::Global_Malloc<T>(diagonal_mapping.descr.size, DataHolder::HOST);
	b_linear_host = AuxFuncCPX::Global_Malloc<T>(diagonal_mapping.descr.size, DataHolder::HOST);
	//if we don't set b_linear_host beyond unknowns, there may be some very large numbers
	//then the convergence test will be wrong
#pragma omp parallel for
	for (int i = 0;i < diagonal_mapping.descr.size;i++) b_linear_host[i] = 0;

	//Initialize preconditioner
	pred.Init(diagonal_mapping.xDoF());

	//Initialize CG system
	cg.preconditioner = &pred;
	cg.linear_mapping = &diagonal_mapping;
	cg.Init(max_iter, relative_tolerance);
}

template<class T, int d>
void DiagonalPoissonSolver<T, d>::Solve(const Field<T, d>& cell_b)
{
	PoissonDescriptor<d>& descr = diagonal_mapping.descr;
	AuxFuncCPX::Fill_CPX_Grid(
		physical_grid.grid,
		[&](const VectorDi& cell) ->T {return cell_b(cell);},
		descr.grid,
		b_linear_host
		);
	AuxFuncCPX::Global_Copy_Array(cg.b_dev, b_linear_host, cg.dof, DataHolder::DEVICE, DataHolder::HOST);
	auto fix_to_zero = [=]__device__(T & b_dev_i, T fixed_dev_i) { if (fixed_dev_i) b_dev_i = 0; };
	cwise_mapping_wrapper(cg.b_dev, descr.d_fixed, fix_to_zero, cg.dof);
	cg.Solve_Device(cg.x_dev, cg.b_dev);
	AuxFuncCPX::Global_Copy_Array(x_linear_host, cg.x_dev, cg.dof, DataHolder::HOST, DataHolder::DEVICE);
	cudaDeviceSynchronize();
	//cg.Solve(x_linear_host, b_linear_host);
	simplexSupport::solver2simplex(descr.grid, x_linear_host, x_field_host);
}


template class DiagonalPoissonSolver<Scalar, 2>;
template class DiagonalPoissonSolver<Scalar, 3>;