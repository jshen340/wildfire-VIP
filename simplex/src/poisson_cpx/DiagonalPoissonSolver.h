//////////////////////////////////////////////////////////////////////////
// CPX solver for a diagonal poisson system
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "DiagonalPoissonMapping.h"
#include "AuxFuncCPX.h"
#include "ConjugatedGradient.h"
#include "Field.h"



template<class T, int d>
class DiagonalPoissonSolver {
	Typedef_VectorDii(d);
public:
	MacGrid<d> physical_grid;
	Field<T, d> x_field_host;
	T* x_linear_host = nullptr, * b_linear_host = nullptr;

	DiagonalPoissonMapping<T, d> diagonal_mapping;
	DiagonalPoissonPreconditioner<T, d> pred;
	ConjugatedGradient cg;
	//padding happens in this level
	//mapping and preconditioner will precisely perform
	//the action that their passed parameters imply
	void Init(const MacGrid<d>& mac_grid, const int max_iter = 50, const Scalar relative_tolerance = 1e-5);
	template<class F1, class F2, class F3>
	void Update_System(F1 is_unknown_func, F2 additional_diag_func, F3 vol_func) {
		//F1 map cell coordinate(VectorDi) to bool
		//F2 map cell coordinate to Scalar
		//F3 map cell coordinate to Scalar
		PoissonDescriptor<d>& descr = diagonal_mapping.descr;
		//fill fixed
		AuxFuncCPX::Fill_CPX_Grid<bool, d>(
			physical_grid.grid,
			[=](const VectorDi& cell)->bool {return !is_unknown_func(cell);},
			descr.grid,
			descr.h_fixed
			);
		//fill additional diagonal term
		AuxFuncCPX::Fill_CPX_Grid<Scalar, d>(
			physical_grid.grid,
			additional_diag_func,
			descr.grid,
			diagonal_mapping.additional_diag_host
			);
		//fill poisson system
		AuxFuncCPX::Fill_CPX_Face<Scalar, d>(physical_grid, vol_func, descr.grid, descr.h_vol);
		diagonal_mapping.To_Device();

		//update preconditioner
		pred.Compute(diagonal_mapping);
	}
	void Solve(const Field<T, d>& cell_b);
};

