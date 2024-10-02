//#####################################################################
// Material point method
// Copyright (c) (2018-), Jiayin Hu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __MaterialPoint_h__
#define __MaterialPoint_h__
#include "Grid.h"
#include "BoundaryCondition.h"
#include "KrylovSolver.h"
#include "Hashtable.h"
#include "GeometricMultiGrid.h"
#include "Particles.h"
#include "Field.h"

template<int d> class MaterialPoint
{
	Typedef_VectorDii(d); Typedef_MatrixD(d);
public:
	const real E = (real)5e1;	////was 1e4
	const real nu = (real)0.2;
	real mu;
	real lambda;

	real vol = (real)1;
	real one_over_dx;
	real one_over_dx_sq;

	Grid<d> grid;
	Particles<d> particles;
	Field<VectorD,d> velocity;
	Field<real,d> mass;
	Field<VectorD,d> momentum;
	Field<VectorD,d> forces;
	Array<MatrixD> deformation_gradient;
	Array<MatrixD> C;

	VectorD g = VectorD::Unit(1)*(real)-9.8;
	bool use_body_force=false;

	virtual void Initialize()
	{
		Initialize_Parameters();
		Initialize_Fields();
	}

	virtual void Initialize_Parameters()
	{
		mu = E / (2 * (1 + nu));
		lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
		one_over_dx=(real)1/grid.dx;
		one_over_dx_sq=pow(one_over_dx,2);
	}
	
	virtual void Initialize_Fields()
	{
		velocity.Resize(grid.cell_counts, VectorD::Zero());
		mass.Resize(grid.cell_counts, (real)1);
		momentum.Resize(grid.cell_counts, VectorD::Zero());
		forces.Resize(grid.cell_counts, VectorD::Zero());
		deformation_gradient.resize(particles.Size(),MatrixD::Identity());
		C.resize(particles.Size(),MatrixD::Zero());
	}

	virtual void Advance(const real dt)
	{
		Particle_To_Grid();
		Update_Grid_Velocities(dt);
		Grid_To_Particle();
		Advection(dt);
		Update_Particle_Deformation_Gradient(dt);
	}

	virtual void Advection(const real dt)
	{
		for (int p = 0; p < particles.Size(); p++) {
			particles.X(p) += particles.V(p) * dt;}
	}
	
	void Particle_To_Grid() 
	{
		mass.Fill(0);
		momentum.Fill(VectorD::Zero());
		////interpolate mass and momentum from particles to grid
		for (int p = 0; p < particles.Size(); p++) {
			const VectorDi cell=grid.Cell_Coord(particles.X(p));
			for(int i=0;i<Grid<d>::Number_Of_Nb_R();i++){
				const VectorDi nb=grid.Nb_R(cell,i);if(!grid.Valid_Cell(nb))continue;
				mass(nb) += particles.M(p) * interpolation_function(grid.Center(nb), particles.X(p));
				momentum(nb) += interpolation_function(grid.Center(nb), particles.X(p)) * particles.M(p) * (particles.V(p) + C[p] * (grid.Center(nb) - particles.X(p)));}}
		
		velocity.Fill(VectorD::Zero());
		iterate_cell(iter, grid) {const VectorDi cell = iter.Coord();
			if (mass(cell) > (real)0) velocity(cell) = momentum(cell) / mass(cell);}
	}

	void Grid_To_Particle() 
	{
		for (int p = 0; p < particles.Size(); p++) {
			C[p] = MatrixD::Zero();
			particles.V(p) = VectorD::Zero();
			const VectorDi cell=grid.Cell_Coord(particles.X(p));
			////interpolate velocity and C from grid to particles
			for(int i=0;i<Grid<d>::Number_Of_Nb_R();i++){
				const VectorDi nb=grid.Nb_R(cell,i);if(!grid.Valid_Cell(nb))continue;
				particles.V(p) += velocity(nb) * interpolation_function(grid.Center(nb), particles.X(p));
				C[p] += velocity(nb) * (grid.Center(nb) - particles.X(p)).transpose() * 
					interpolation_function(grid.Center(nb), particles.X(p)) * 4*one_over_dx_sq;}}
	}

	void Update_Grid_Velocities(const real dt) 
	{
		forces.Fill(VectorD::Zero());
		////update force using MLS-MPM
		for (int p = 0; p < particles.Size(); p++) {
			const VectorDi cell=grid.Cell_Coord(particles.X(p));
			////update grid force from particle F
			for(int i=0;i<Grid<d>::Number_Of_Nb_R();i++){
				const VectorDi nb=grid.Nb_R(cell,i);if(!grid.Valid_Cell(nb))continue;
				forces(nb) -= interpolation_function(grid.Center(nb), particles.X(p)) * vol * 4 * one_over_dx_sq * 
					PK1_Neohookean(deformation_gradient[p]) * deformation_gradient[p].transpose() * (grid.Center(nb) - particles.X(p));}}
		////time integration
		for (int i = 0; i < grid.Number_Of_Cells(); i++) {
			forces(i) += mass(i) * g;
			if (mass(i) > 0) velocity(i) += dt * forces(i) / mass(i);
			for (int axis = 0; axis < d; axis++){
				if (grid.Cell_Coord(i)[axis] ==0 || grid.Cell_Coord(i)[axis]==grid.cell_counts[axis] - 1)
					velocity(i)[axis] = (real)0;}}
	}

	void Update_Particle_Deformation_Gradient(const real dt) 
	{
		for (int p = 0; p < particles.Size(); p++) {
			deformation_gradient[p] = (MatrixD::Identity() + dt * C[p]) * deformation_gradient[p];}
	}

protected:
	real kernel(real x1, real x2) 
	{ //quadratic kernel
		real x = abs(x2 - x1);
		if (0 <= x && x < 0.5) return 0.75 - x * x;
		else if (0.5 <= x && x < 1.5) return 0.5 * (1.5 - x) * (1.5 - x);
		else return 0;
	}

	real kernel_gradient(real x1, real x2) 
	{
		real x = x2 - x1;
		if (-0.5 < x && x < 0.5) return -2 * x;
		else if (-1.5 < x && x <= -0.5) return 1.5 + x;
		else if (0.5 <= x && x < 1.5) return -1.5 + x;
		else return 0;
	}

	real interpolation_function(VectorD p1, VectorD p2) 
	{
		real s = 1;
		for (int i = 0; i < d; i++) {
			s *= kernel(p1[i]*one_over_dx, p2[i]*one_over_dx);}
		return s;
	}

	VectorD interpolation_function_gradient(VectorD p1, VectorD p2) 
	{
		VectorD s = VectorD::Ones();
		for (int i = 0; i < d; i++) {
			for (int j = 0; j < d; j++) {
				if (j == i) s[i] *= kernel_gradient(p1[i]*one_over_dx, p2[i]*one_over_dx)*one_over_dx;
				else s[i] *= kernel(p1[j]*one_over_dx, p2[j]*one_over_dx);}}
		return s;
	}

	MatrixD PK1_Neohookean(const MatrixD& F)
	{
		real J = F.determinant(); MatrixD F_invT = F.inverse().transpose(); 
		auto pk1 = mu * (F - F_invT) + lambda * log(J) * F_invT;
		return pk1;
	}

	MatrixD PK1_Corotated(const MatrixD& F)
	{
		real J = F.determinant(); MatrixD F_invT = F.inverse().transpose();
		Eigen::JacobiSVD<MatrixD> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);
		MatrixD R = svd.matrixU() * svd.matrixV().transpose();
		auto pk1 = 2 * mu * (F - R) + lambda * (J - 1) * J * F_invT;
		return pk1;
	}
};
#endif