//////////////////////////////////////////////////////////////////////////
// Use sparse solver to solve Poisson systems
// The reason is that some Poisson systems has very few unknowns against a large computational field
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "SparseSolverCPX.h"
#include "Field.h"
#include "MacGrid.h"
#include "ArrayFunc.h"

template<class T, int d>
class SparsePoissonSolverCPX : public SparseSolverCPX<T> {
	using Base = SparseSolverCPX<T>;
	using Base::sparse_mapping;
	using Base::pred;
	using Base::cg;
public:
	Typedef_VectorDii(d);
	//for poisson system
	int dof;
	MacGrid<d> mac_grid;
	Array<int> c_idx;
	Field<int, d> cell_idx;
	Array<VectorDi> idx_coord;
	Array<int> idx_row_nnz;
	Array<real> rhs;
	Array<real> x_linear;
	Field<real, d> x_field;
	void Init(const MacGrid<d>& _mac_grid, std::function<bool(const VectorDi&)> is_unknown_func, std::function<Scalar(const int, const VectorDi&)> vol_func, const int max_iter = 50, const Scalar relative_tolerance = 1e-5){
		mac_grid = _mac_grid;
		//we call the index in serialized matrix "idx"
		//the coordinate of face node is "face_cell"
		//the serial id in face_grid of a face node is "c"
		const int cell_num = mac_grid.grid.Number_Of_Cells();
		//temporary structure for boosting, c_idx[c]=1 or corresponding idx
		c_idx.resize(cell_num);
		mac_grid.grid.Exec_Each(
			[&](const VectorDi& cell) {
				int c = mac_grid.grid.Cell_Index(cell);
				if (is_unknown_func(cell)) c_idx[c] = 1;
				else c_idx[c] = -1;
			}
		);
		//Indexing every cell, getting number of valid elements
		dof = 0;
		for (int i = 0;i < cell_num;i++) {
			if (c_idx[i] != -1) {
				c_idx[i] = dof++;
			}
		}
		//data structure for cell-idx mapping
		cell_idx.Resize(mac_grid.grid.cell_counts);
		idx_coord.resize(dof);
		mac_grid.grid.Exec_Each(
			[&](const VectorDi& cell) {
				int c = mac_grid.grid.Cell_Index(cell);
				int idx = c_idx[c];
				if (idx != -1) {
					cell_idx(cell) = idx;
					idx_coord[idx] = cell;
				}
				else cell_idx(cell) = -1;
			}
		);
		//calculate non-zero numbers, each line and total
		idx_row_nnz.resize(dof);
		ArrayFunc::Calc_Each(
			[&](const int idx) {
				int row_nnz = 1;//diagonal element itself
				const VectorDi cell = idx_coord[idx];
				for (int i = 0;i < Grid<d>::Number_Of_Nb_C();i++) {
					VectorDi nb_cell = Grid<d>::Nb_C(cell, i);
					if (mac_grid.grid.Valid_Cell(nb_cell) && cell_idx(nb_cell) != -1) {
						row_nnz++;
					}
				}
				return row_nnz;
			},
			idx_row_nnz
				);
		int non_zero_numbers = 0;
#pragma omp parallel for reduction(+:non_zero_numbers)
		for (int idx = 0;idx < dof;idx++) {
			non_zero_numbers += idx_row_nnz[idx];
		}
		Eigen::SparseMatrix<T, Eigen::RowMajor, int> A;
		
		A.resize(dof, dof);
		A.resizeNonZeros(non_zero_numbers);

		//[SERIAL] iterate all elements, calculate a prefix sum
		int* outer = A.outerIndexPtr();outer[0] = 0;
		for (int idx = 0;idx < dof;idx++) {
			outer[idx + 1] = outer[idx] + idx_row_nnz[idx];
		}
		ArrayFunc::Exec_Each(
			[&](const int idx) {
				VectorDi cell = idx_coord[idx];
				T dia_coef = 0;
				Array<std::pair<int, T>> col_val_pairs;col_val_pairs.clear();
				for (int i = 0;i < MacGrid<d>::Number_Of_Cell_Incident_Faces();i++) {
					int axis;VectorDi face;
					MacGrid<d>::Cell_Incident_Face(cell, i, axis, face);
					VectorDi nb_cell = Grid<d>::Nb_C(cell, i);
					T a = vol_func(axis, face);
					dia_coef += a;
					if (mac_grid.grid.Valid_Cell(nb_cell)) {
						int col = cell_idx(nb_cell);
						if (col != -1) {
							col_val_pairs.push_back(std::make_pair(col, -a));
						}
					}
				}
				col_val_pairs.push_back(std::make_pair(idx, dia_coef));
				int* col = A.innerIndexPtr() + outer[idx];
				T* val = A.valuePtr() + outer[idx];
				std::sort(col_val_pairs.begin(), col_val_pairs.end());
				for (int k = 0;k < col_val_pairs.size();k++) {
					col[k] = col_val_pairs[k].first;
					val[k] = col_val_pairs[k].second;
				}
			},
			idx_coord
				);

		sparse_mapping.Update(A);
		pred.Init(sparse_mapping);
		cg = ConjugatedGradient();
		cg.preconditioner = &pred;
		cg.linear_mapping = &sparse_mapping;
		cg.Init(max_iter, relative_tolerance);
		checkCudaErrors(cudaGetLastError());
	}
	void Fill_RHS(std::function<T(const VectorDi&)> rhs_func) {
		rhs.resize(dof);
		ArrayFunc::Calc_Each(
			[&](const int idx) {
				return rhs_func(idx_coord[idx]);
			},
			rhs
				);
	}
	void Solve(void) {
		x_linear.resize(dof);
		Base::Solve(x_linear.data(), rhs.data());
		x_field.Resize(mac_grid.grid.cell_counts, (T)0);
		mac_grid.grid.Exec_Each(
			[&](const VectorDi& cell) {
				if (cell_idx(cell) != -1) {
					x_field(cell) = x_linear[cell_idx(cell)];
				}
				else x_field(cell) = 0;
			}
		);
	}
};

