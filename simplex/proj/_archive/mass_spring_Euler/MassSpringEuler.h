//#####################################################################
// mass spring Euler
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __MassSpringEuler_h__
#define __MassSpringEuler_h__
#include "Common.h"
#include "Hashtable.h"
#include "Particles.h"
#include "Grid.h"
#include "Field.h"
#include "Interpolation.h"

template<int d> class MassSpringEuler
{
	Typedef_VectorDii(d); Typedef_MatrixD(d);
public:
	real ks_0 = (real)1;
	real kd_0 = (real)1;
	real rest_length=(real)0.25;

	Particles<d> particles;
	Grid<d> grid;
	Field<real,d> potential;

	Hashtable<int, VectorD> psi_D_values;	////displacements
	Hashtable<int, VectorD> psi_N_values;	////forces
	VectorD g = VectorD::Unit(1) * (real)-1.;
	bool use_body_force = true;

	virtual void Initialize() {
		potential.Resize(grid.cell_counts);
	}

	virtual void Advance(const real dt)
	{
		Update_Field();
		Update_Forces();
		Enforce_Boundary_Conditions(particles.VRef(), particles.FRef());
		for (int i = 0; i < particles.Size(); i++) {
			particles.V(i) *= (real).95;
			particles.V(i) += particles.F(i) / particles.M(i) * dt;
			//particles.V(i) *= (real).95;
			particles.X(i) += particles.V(i) * dt;
		}
	}

	void Update_Field() {
		potential.Fill((real)0.);
		for (int p = 0; p < particles.Size(); p++) {
			for (int i = 0; i < potential.array.size(); i++) {
				//potential(i) += (real).5 * ks_0 * pow((grid.Center(grid.Cell_Coord(i)) - particles.X(p)).norm() - rest_length, 2);
				//potential(i) += (real).5 * ks_0 * (grid.Center(grid.Cell_Coord(i)) - particles.X(p)).squaredNorm();
				//potential(i) += -(real)0.01 / (grid.Center(grid.Cell_Coord(i)) - particles.X(p)).norm();
				//potential(i) += pow((grid.Center(grid.Cell_Coord(i)) - particles.X(p)).norm() - rest_length, 4);
				potential(i) += pow(abs(grid.Center(grid.Cell_Coord(i))[0]-particles.X(p)[0])-rest_length,2) + pow(abs(grid.Center(grid.Cell_Coord(i))[1]-particles.X(p)[1])-rest_length,2);
				//potential(i) += pow(abs(grid.Center(grid.Cell_Coord(i))[0] - particles.X(p)[0]) - rest_length, 2);
				//potential(i) += sin((grid.Center(grid.Cell_Coord(i)) - particles.X(p)).norm()-rest_length);
				//potential(i) += (grid.Center(grid.Cell_Coord(i)) - particles.X(p))[0] * (grid.Center(grid.Cell_Coord(i)) - particles.X(p))[1];
				//potential(i) += (real).5 * ks_0 * pow((grid.Center(grid.Cell_Coord(i)) - particles.X(p)).norm() - rest_length, 2) + pow(abs(grid.Center(grid.Cell_Coord(i))[0] - particles.X(p)[0]) - rest_length, 2);
				//potential.array[i] += log(1 + pow((grid.Center(grid.Cell_Coord(i)) - particles.X(p)).norm() - rest_length, 2));
			}
		}
	}

	void Update_Forces()
	{
		Field<VectorD, d> gradients(grid.cell_counts, VectorD::Zero());
		iterate_cell(iter, grid) {
			const VectorDi cell = iter.Coord();
			for (int axis = 0; axis < d; axis++) {
				VectorDi nb_left = cell - VectorDi::Unit(axis);
				if (!grid.Valid_Cell(nb_left))nb_left = cell;
				VectorDi nb_right = cell + VectorDi::Unit(axis);
				if (!grid.Valid_Cell(nb_right))nb_right = cell;
				gradients(cell)[axis] = (potential(nb_right) - potential(nb_left)) / ((real)2* grid.dx);
			}
		}

		Interpolation<d> intp(grid);
		for (int i = 0; i < particles.Size(); i++) {
			VectorD grad=intp.Interpolate_Centers(gradients, particles.X(i));
			particles.F(i) = -grad;
		}

		if (use_body_force)
			for (int i = 0; i < particles.Size(); i++) {
				particles.F(i) += particles.M(i) * g;
			}
		for (auto p : psi_N_values) { int idx = p.first; const VectorD& f = p.second; particles.F(idx) += f; }
	}

	void Set_Psi_D(const int p, const VectorD v = VectorD::Zero()) { psi_D_values[p] = v; }
	bool Is_Psi_D(const int p) { return psi_D_values.find(p) != psi_D_values.end(); }
	void Set_Psi_N(const int p, const VectorD f) { psi_N_values[p] = f; }
	bool Is_Psi_N(const int p) { return psi_N_values.find(p) != psi_N_values.end(); }

	void Enforce_Boundary_Conditions(Array<VectorD>& V, Array<VectorD>& F)
	{
		for (auto p : psi_D_values) { int idx = p.first; const VectorD& v = p.second; V[idx] = v; F[idx] = VectorD::Zero(); }
	}
};

#endif