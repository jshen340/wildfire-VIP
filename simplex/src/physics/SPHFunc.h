//////////////////////////////////////////////////////////////////////////
// Basic algorithms for different SPH methods
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Common.h"
#include "Kernels.h"
#include "Particles.h"
#include "SPHParticles.h"
#include "Integrators.h"

namespace SPHFunc {

	//VS for V times Scalar
	template<typename T1, typename T2>
	decltype(auto) MultiplyVS_Func(const Array<T1> &arr, const T2 &alpha){
		return [&](const int& i) {
			return arr[i] * alpha;
		};
	}
	template<class T>
	SingleIFunc<T> Index_Function(const Array<T>& arr) {
		return [&](const int& idx) {return arr[idx]; };
	}

	//mass density = m_i * number density
	template<int d>
	real Rest_Number_Density(const real& dx, const KernelSPH& kernel);

	//template<int d>
	//int IISPH_Projection(const Array<T> &)

	//the returned coeff alpha is used in PCI algorithm.
	//say, predicted pressure p_i=beta*(rho_0-rho_i), where beta=alpha/dt^2. 
	template<int d>	
	real Calculate_PCI_Coefficient(const real dx, const KernelSPH &kernel);

	//return a f(i,j), where i is the point we're taking care of
	template<int d>	PairIFuncV<d> Symmetric_Grad_Term_Func(const SPHParticles<d>& particles, const Array<real>& arr, const KernelSPH& kernel);
	template<int d>	SingleIFuncV<d> Symmetric_Grad_Func(const SPHParticles<d>& particles, const Array<real>& arr, const KernelSPH& kernel);

	template<int d>	PairIFuncV<d> Viscosity_Laplacian_Term_Func(const SPHParticles<d>& particles, const Array<Vector<real, d> >& arr, const KernelSPH& kernel);
	template<int d>	SingleIFuncV<d> Viscosity_Laplacian_Func(const SPHParticles<d>& particles, const Array<Vector<real, d> >& arr, const KernelSPH& kernel);

	template<int d>	PairIFuncS Mass_Density_Term_Func(const Array<Vector<real, d> >& x, const Array<real>& m, const KernelSPH& kernel);
	template<int d>	SingleIFuncS Mass_Density_Func(const SPHParticles<d>& particles, const KernelSPH& kernel);

	//p = kp * ((rho / rho_0) ^ kb - 1)
	template<int d>
	SingleIFuncS State_Equation(const SPHParticles<d>& particles, const real& rest_density, const real& kp, const int& kb);
}

//Predictive-Corrective-Incompressible SPH
//The parameter "mode" may be different from that in particles. Because the former is for gradient calculation, and the latter is for mass density calculation.
template<int d>
class PCISPHSolver {
	Typedef_VectorDii(d);
public:
	std::shared_ptr<IntegratorX<d>> intg_x_ptr;
	real rho_0;
	real alpha;
	KernelSPH grad_kernel;//seems that only SPIKY kernel works
	int max_iter;
	real max_rel_err;
	bool verbose = false;
	PCISPHSolver(std::shared_ptr<IntegratorX<d>> _intg, real _rho_0, real _alpha, const KernelSPH &_kernel, int _max_iter = 20, real _max_rel_err = 1e-3) {
		intg_x_ptr = _intg;
		rho_0 = _rho_0;
		alpha = _alpha;
		grad_kernel = _kernel;
		max_iter = _max_iter;
		max_rel_err = _max_rel_err;
	}
	//just modify velocity. does not perform time integration
	real Solve(SPHParticles<d>& particles, const real& dt);
};

template<int d>
class IISPHSolver {
	Typedef_VectorDii(d);
public:
	std::shared_ptr<IntegratorX<d>> intg_x_ptr;
	real rho_0;
	real omega;//realxation coefficient
	int max_iter;
	real max_rel_err;
	bool verbose = false;
	IISPHSolver(std::shared_ptr<IntegratorX<d>> _intg, real _rho_0, real _omega, int _max_iter = 20, real _max_rel_err = 1e-3) {
		intg_x_ptr = _intg;
		rho_0 = _rho_0;
		omega = _omega;
		max_iter = _max_iter;
		max_rel_err = _max_rel_err;
	}
	//just change velocity with respect to pressure
	real Solve(SPHParticles<d>& particles, const real& dt);
};


