//////////////////////////////////////////////////////////////////////////
// Poisson solver interface aligning with simplex for CPX solver
// Copyright (c) (2018-), Yueyang Xianzang, Jinyuan Liu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#pragma once
#include "Grid.h"
#include "Field.h"
#include "FaceField.h"
#include "PoissonDescriptor.h"
#include "ConjugatedGradient.h"

template<int d>
class Poisson
{
	Typedef_VectorDii(d);
public:
	const int block_size = (d == 2 ? 8 : 4);
	//Grid<d> grid;
	MacGrid<d> physical_grid;//original, unpadded grid from simplex solver
	int l=5;

	//mg_descr describes A for equation Ax=b
	PoissonDescriptor<d> *mg_descr;
public:
	ConjugatedGradient cg;
	bool closed = false;

	Field<Scalar, d> x;//, b;
	Scalar *temp_x, *temp_b;

public:
	Poisson() {
		//if (d == 2) {
		//	grid_size = VectorDi::Ones() * 128;
		//	l = 5;
		//}
		//else {
		//	grid_size = VectorDi::Ones() * 32;
		//	l = 4;
		//}
	}

	//default value for Eigen is (dof*2) iterations and std::numeric_limits<T>:: epsilon()
	//void init(VectorDi _grid_size, Scalar _dx, const int max_iter = 50, const Scalar relative_tolerance = 1e-5);
	void Init(const MacGrid<d>& mac_grid, int max_iter = -1, Scalar relative_tolerance = -1);

	void Solve();

	void Solve_Fast();

	void Init_Boundary(const FaceField<Scalar, d>& face_vol, const Field<int, d>& cell_fixed, bool _closed);

	void Update_b(const Field<Scalar, d>& cell_b);

	//is_unknown_func(cell)==true means this cell is a unknown variable in system, and its fixed is set to false
	void Update_Unknown(std::function<bool(const VectorDi&)> is_unknown_func);

	void Update_Vol(std::function<Scalar(const int, const VectorDi&)> vol_func);

	void Send_To_Device(void);

	template<class Fcell> //Fcell: T(const VectorDi &cell)
	void Update_RHS(Fcell rhs_func) {
		PoissonDescriptor<d>& descr = mg_descr[l - 1];
		int cell_num = physical_grid.grid.Number_Of_Cells();
#pragma omp parallel for
		for (int i = 0; i < cell_num; i++) {
			VectorDi cell = physical_grid.grid.Cell_Coord(i);
			int cell_ind = -1;
			if constexpr (d == 2) cell_ind = descr.grid.cell_ind(cell[0], cell[1]);
			else if constexpr (d == 3) cell_ind = descr.grid.cell_ind(cell[0], cell[1], cell[2]);
			temp_b[cell_ind] = rhs_func(cell);
		}
	}
};
