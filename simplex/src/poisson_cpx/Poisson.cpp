#include "Poisson.h"

#include "PoissonMapping.h"
#include "PoissonMapping3D.h"
#include "PoissonLike.h"
#include "PoissonSubSystem.h"
#include "PoissonSubSystem3D.h"
#include "MultigridSolve.h"

#include "SimplexSupport.h"

#include "Timer.h"

template<int d>
void Poisson<d>::Init(const MacGrid<d>& mac_grid, int max_iter, Scalar relative_tolerance)
{
	physical_grid = mac_grid;

	//padding to multiple of 8x8 in 2d, 4x4x4 in 3d
	VectorDi grid_size = AuxFuncCPX::Round_Up_To_Align<d>(physical_grid.grid.cell_counts, block_size);
	int n = grid_size.minCoeff();
	if (d == 2) l = round(std::log2(n & -n)) - 2;
	else l = round(std::log2(n & -n)) - 1;

	x.Resize(grid_size);
	//b.Resize(grid_size);

	PoissonDescriptor<d> descr = PoissonDescriptor<d>();
	descr.init(grid_size);
#pragma omp parallel for
	for (int i = 0;i < descr.fsize;i++) descr.h_vol[i] = (Scalar)0;
#pragma omp parallel for
	for (int i = 0;i < descr.size;i++) descr.h_fixed[i] = true;

	temp_x = new Scalar[descr.size];
	temp_b = new Scalar[descr.size];

	//if we don't set temp_b beyond unknowns, there may be some very large numbers
	//then the convergence test will be wrong
#pragma omp parallel for
	for (int i = 0;i < descr.size;i++) temp_b[i] = 0;

	MultigridSolve *mg_solve = new MultigridSolve();
	mg_solve->init(l);

	mg_descr = new PoissonDescriptor<d>[l];
	mg_descr[l - 1] = descr;

	for (int i = l - 2; i >= 0; i--) mg_descr[i] = createSubSystem(mg_descr[i + 1]);

	using FixedMapping = typename If<d == 2, PoissonMappingFixed, PoissonMapping3DFixed >::Type;
	using DownSampler= typename If<d == 2, PoissonDownSample, PoissonDownSample3D >::Type;
	using UpSampler= typename If<d == 2, PoissonUpSample, PoissonUpSample3D >::Type;
	for (int i = 0; i < l; i++)
	{
		FixedMapping*t0 = new FixedMapping();
		t0->init(&mg_descr[i]);
		mg_solve->mapping[i] = t0;
		if (i != 0)
		{
			//PoissonSmoothing *t1 = new PoissonSmoothing(t0);
			LinearMapping *t1 = PoissonLike::GenerateJacobiPreconditioner(mg_descr[i].grid, t0);
			mg_solve->preSmoother[i] = t1;
			//mg_solve->preSmoother[i] = nullptr;

			//PoissonSmoothing *t2 = new PoissonSmoothing(t0);
			mg_solve->postSmoother[i] = t1;
			//mg_solve->postSmoother[i] = nullptr;

			DownSampler*t3 = new DownSampler(&mg_descr[i], &mg_descr[i - 1]);
			mg_solve->downSample[i - 1] = t3;

			UpSampler*t4 = new UpSampler(&mg_descr[i], &mg_descr[i - 1]);
			mg_solve->upSample[i - 1] = t4;
		}
		else
		{
			//PoissonDirectSolve *direct = new PoissonDirectSolve(t0);
			//direct->finish();
			LinearMapping *direct = PoissonLike::GenerateDirectSolve(mg_descr[i].grid, t0, closed); //default value is closed=false
			//LinearMapping *direct = PoissonLike::GenerateDirectSolve(mg_descr[i].grid, t0, false); //was closed originally for false
			mg_solve->preSmoother[i] = direct;
			mg_solve->postSmoother[i] = nullptr;
		}
	}

	mg_solve->finish();

	cg = ConjugatedGradient();
	//cg.relative_tolerance = 1e-5;
	//cg.linear_mapping = new PoissonDiagMapping(mapping);
	//cg.preconditioner = new DescentSmoother(new PoissonSmoothing((PoissonMapping*)cg.linear_mapping), 20, cg.linear_mapping);
	cg.preconditioner = mg_solve;
	cg.linear_mapping = mg_solve->mapping[l - 1];
	if (max_iter == -1) { max_iter = 2*descr.size; } //default value
	if (relative_tolerance == (Scalar) -1) { relative_tolerance=std::numeric_limits<Scalar>::epsilon(); } //default value
	cg.Init(max_iter, relative_tolerance);
}

template<int d>
void Poisson<d>::Solve()
{
	for (int i = l - 2; i >= 0; i--) updateSubSystem(mg_descr[i], mg_descr[i + 1]);
	cg.preconditioner->update();

	cudaDeviceSynchronize();
	checkCudaErrors(cudaGetLastError());

	PoissonDescriptor<d>& descr = mg_descr[l - 1];

	//simplexSupport::simplex2solver(b, descr.grid, temp_b);

	cg.Solve(temp_x, temp_b);

	simplexSupport::solver2simplex(descr.grid, temp_x, x);
}

template<int d>
void Poisson<d>::Solve_Fast()
{
	for (int i = l - 2; i >= 0; i--) updateSubSystem(mg_descr[i], mg_descr[i + 1]);
	cg.preconditioner->update();

	cudaDeviceSynchronize();
	checkCudaErrors(cudaGetLastError());

	PoissonDescriptor<d>& descr = mg_descr[l - 1];

	//simplexSupport::simplex2solver(b, descr.grid, temp_b);

	cg.Solve(temp_x, temp_b);

	int cell_num = physical_grid.grid.Number_Of_Cells();
#pragma omp parallel for
	for (int i = 0; i < cell_num; i++) {
		VectorDi cell = physical_grid.grid.Cell_Coord(i);
		int cell_ind = -1;
		if constexpr (d == 2) cell_ind = descr.grid.cell_ind(cell[0], cell[1]);
		else if constexpr (d == 3) cell_ind = descr.grid.cell_ind(cell[0], cell[1], cell[2]);
		x(cell) = temp_x[cell_ind];
	}
	//simplexSupport::solver2simplex(descr.grid, temp_x, x);
}

template<int d>
void Poisson<d>::Init_Boundary(const FaceField<Scalar, d>& face_vol, const Field<int, d>& cell_fixed, bool _closed)
{
	PoissonDescriptor<d>& descr = mg_descr[l - 1];

	int size = descr.size;
	bool *fixed = new bool[size];
	Scalar* vol = new Scalar[descr.fsize];
	memset(fixed, 0, sizeof(bool) * descr.size);
	memset(vol, 0, sizeof(Scalar) * descr.fsize);

	simplexSupport::simplex2solver(cell_fixed, descr.grid, fixed, [=](int v)->bool { if (v) return true; else return false; });

	simplexSupport::simplex2solver(face_vol, descr.grid, vol);

	descr.setFixed(fixed);
	descr.setVol(vol);
	descr.toDevice();
	descr.finish();

	closed = _closed;

	delete[] fixed;
	delete[] vol;
}

template<int d>
void Poisson<d>::Update_b(const Field<Scalar, d>& cell_b)
{
	PoissonDescriptor<d>& descr = mg_descr[l - 1];
	AuxFuncCPX::Fill_CPX_Grid<Scalar, d>(
		physical_grid.grid,
		[&](const VectorDi& cell) {return cell_b(cell);},
		descr.grid,
		temp_b
	);
//	int cell_num = physical_grid.grid.Number_Of_Cells();
//#pragma omp parallel for
//	for (int i = 0; i < cell_num; i++) {
//		VectorDi cell = physical_grid.grid.Cell_Coord(i);
//		int cell_ind = -1;
//		if constexpr (d == 2) cell_ind = descr.grid.cell_ind(cell[0], cell[1]);
//		else if constexpr (d == 3) cell_ind = descr.grid.cell_ind(cell[0], cell[1], cell[2]);
//		temp_b[cell_ind] = cell_b(cell);
//	}
	//b = cell_b;
}

template<int d>
void Poisson<d>::Update_Unknown(std::function<bool(const VectorDi&)> is_unknown_func)
{
	PoissonDescriptor<d>& descr = mg_descr[l - 1];
	AuxFuncCPX::Fill_CPX_Grid<bool, d>(
		physical_grid.grid,
		[=](const VectorDi& cell)->bool {return !is_unknown_func(cell);},
		descr.grid,
		descr.h_fixed
	);
//	PoissonDescriptor<d>& descr = mg_descr[l - 1];
//	int cell_num = physical_grid.grid.Number_Of_Cells();
//#pragma omp parallel for
//	for (int i = 0; i < cell_num; i++) {
//		VectorDi cell = physical_grid.grid.Cell_Coord(i);
//		int cell_ind = -1;
//		if constexpr (d == 2) cell_ind = descr.grid.cell_ind(cell[0], cell[1]);
//		else if constexpr (d == 3) cell_ind = descr.grid.cell_ind(cell[0], cell[1], cell[2]);
//		descr.h_fixed[cell_ind] = !is_unknown_func(cell);
//	}
}

template<int d>
void Poisson<d>::Update_Vol(std::function<Scalar(const int, const VectorDi&)> vol_func)
{
	PoissonDescriptor<d>& descr = mg_descr[l - 1];
	AuxFuncCPX::Fill_CPX_Face(physical_grid, vol_func, descr.grid, descr.h_vol);
//	int offset = 0;
//	for (int axis = 0; axis < d; axis++) {
//		int face_num = physical_grid.Number_Of_Faces(axis);
//#pragma omp parallel for
//		for (int i = 0; i < face_num; i++) {
//			VectorDi face = physical_grid.Face_Coord(axis, i);
//			int face_ind = -1;
//			if constexpr (d == 2) { face_ind = offset + descr.grid.face_ind(face[0], face[1], axis); }
//			else if constexpr (d == 3) { face_ind = offset + descr.grid.face_ind(face[0], face[1], face[2], axis); }
//			descr.h_vol[face_ind] = vol_func(axis, face);
//		}
//		offset += descr.grid.face_size(axis);
//	}
}

template<int d>
void Poisson<d>::Send_To_Device(void)
{
	PoissonDescriptor<d>& descr = mg_descr[l - 1];
	descr.toDevice();
	descr.finish();
}


template class Poisson<2>;
template class Poisson<3>;