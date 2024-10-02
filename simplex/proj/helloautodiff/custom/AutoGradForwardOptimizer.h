#ifndef __AutoGradForwardOptimizer_h__
#define __AutoGradForwardOptimizer_h__
#include "OptimizerMma.h"
#include "OptimizerIpOpt.h"
#include <iostream>
// Eigen includes
#include <Eigen/Core>

// autodiff include
#include <autodiff/forward.hpp>
#include <autodiff/forward/eigen.hpp>
using namespace autodiff;
#define USE_MMA_FOR_MESH_OPT 1

#ifdef USE_MMA_FOR_MESH_OPT
class AutoGradForwardOptimizer : OptimizerMMA
#else
class AutoGradForwardOptimizer : OptimizerIpOpt
#endif
{
public:
	const int N_VAR = 10;
	VectorXdual x;

	//////////////////////////////////////////////////////////////////////////
	////MMA Modified functions
	void Initialize_Optimizer()
	{
		n_var = N_VAR;
		n_cons = 0;

		var_lb = (real).0;
		var_ub = (real)1.;
		Allocate_Data();

		////parameters for MMA only
		#ifdef USE_MMA_FOR_MESH_OPT
			movlim = (real).1;
			tol = (real)1e-4;
		#endif

		x.resize(N_VAR);
		for (int i = 0; i < N_VAR; i++) {
			x[i] = (real) i*0.1;
		}
		std::cout << x.transpose() << std::endl;

		for (int i = 0; i < N_VAR; i++) {
			var[i] = x[i].val;
		}
	}

	// function has to be static
	static dual Func(const VectorXdual& x)
	{
		return x.cwiseProduct(x).sum(); // sum([x(i) * x(i) for i = 1:10])
	}

	virtual real Compute_Objective(const real* var)
	{
		for (int i = 0; i < N_VAR; i++) {
			x[i] = var[i];
		}

		real obj = Func(x).val;
		std::cout << "current x is: "<<x.transpose() << std::endl;
		std::cout << "current obj is: "<<obj << std::endl;
		return obj;
	}

	virtual void Compute_Gradient(const real* var, real* grad)
	{
		dual u;  // the output scalar u = f(x) evaluated together with gradient below
		VectorXd g = gradient(Func, wrt(x), at(x), u);  // evaluate the function value u and its gradient vector g = du/dx
		for (int i = 0; i < N_VAR; i++) {
			grad[i] = g[i];
		}
	}

	virtual void Optimize()
	{
	#ifdef USE_MMA_FOR_MESH_OPT
		OptimizerMMA::Optimize();
	#else
		OptimizerIpOpt::Optimize();
	#endif
	}

	virtual void Write_Substep(const int frame)
	{
		std::cout << "Optimization iteration " << frame << std::endl;
	}
};
#endif