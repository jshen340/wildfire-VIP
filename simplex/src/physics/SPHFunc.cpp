#include "SPHFunc.h"
#include "PointSetFunc.h"
#include "ArrayFunc.h"

namespace SPHFunc {
	template<int d>
	real Rest_Number_Density(const real& dx, const KernelSPH& kernel)
	{
		Typedef_VectorD(d);
		Points<d> points;
		points.Add_Element();
		points.X(0) = VectorD::Zero();
		if constexpr (d == 2) PointSetFunc::Extend_Hexagon_Packing(points, kernel.h, dx);
		else PointSetFunc::Extend_FCC_Packing(points, kernel.h, dx);
		real rho_0 = 0;
		int N = points.Size();
		for (int i = 0; i < N; i++) rho_0 += kernel.Weight<d>(points.X(0) - points.X(i));
		return rho_0;
	}
	template real Rest_Number_Density<2>(const real&, const KernelSPH&);
	template real Rest_Number_Density<3>(const real&, const KernelSPH&);

	template<int d>
	real Calculate_PCI_Coefficient(const real dx, const KernelSPH& kernel)
	{
		Typedef_VectorD(d);
		Points<d> points;
		points.Add_Element();
		points.X(0) = VectorD::Zero();
		if constexpr (d == 2) PointSetFunc::Extend_Hexagon_Packing(points, kernel.h, dx);
		else PointSetFunc::Extend_FCC_Packing(points, kernel.h, dx);
		real rho_0 = 0;//assume m=1
		VectorD grad_sum = VectorD::Zero();
		real dot_sum = 0;
		int N = points.Size();
		for (int i = 0; i < N; i++) {
			VectorD r = points.X(0) - points.X(i);
			rho_0 += kernel.Weight<d>(r);
			VectorD grad = kernel.Grad<d>(r);
			grad_sum += grad;
			dot_sum += grad.squaredNorm();
		}
		real S = grad_sum.squaredNorm() + dot_sum;
		return -rho_0 * rho_0 / (2 * S);
	}
	template real Calculate_PCI_Coefficient<2>(const real, const KernelSPH&);
	template real Calculate_PCI_Coefficient<3>(const real, const KernelSPH&);

	template<int d>
	PairIFuncV<d> Symmetric_Grad_Term_Func(const SPHParticles<d>& particles, const Array<real>& arr, const KernelSPH& kernel)
	{
		return [&](const int& i, const int& j) {
			const real& rho_i = particles.Rho(i);
			const real& rho_j = particles.Rho(j);
			const real& A_i = arr[i];
			const real& A_j = arr[j];
			Vector<real, d> grad = kernel.Grad<d>(particles.X(i) - particles.X(j));
			return rho_i * particles.M(j) * (A_i / (rho_i * rho_i) + A_j / (rho_j * rho_j)) * grad;
		};
	}

	template<int d>	
	SingleIFuncV<d> Symmetric_Grad_Func(const SPHParticles<d>& particles, const Array<real>& arr, const KernelSPH& kernel) {
		return particles.Neighbor_Sum_Wrapper(Symmetric_Grad_Term_Func(particles, arr, kernel));
	}
	template SingleIFuncV<2> Symmetric_Grad_Func<2>(const SPHParticles<2>&, const Array<real>&, const KernelSPH&);
	template SingleIFuncV<3> Symmetric_Grad_Func<3>(const SPHParticles<3>&, const Array<real>&, const KernelSPH&);

	template<int d>
	PairIFuncV<d> Viscosity_Laplacian_Term_Func(const SPHParticles<d>& particles, const Array<Vector<real, d>>& arr, const KernelSPH& kernel)
	{
		return [&](const int& i, const int& j)->Vector<real, d> {
			if (i == j) return Vector<real, d>::Zero();
			Vector<real, d> rij = particles.X(i) - particles.X(j);
			Vector<real, d> grad = kernel.Grad<d>(rij);
			real lij2 = rij.squaredNorm();
			real one_over_lij2 = 1.0 / (lij2 + 0.001 * kernel.h * kernel.h);
			return  2 * (d + 2) * particles.M(j) / particles.Rho(j) * rij.dot(arr[i] - arr[j]) * one_over_lij2 * grad;
		};
	}
	template<int d>
	SingleIFuncV<d> Viscosity_Laplacian_Func(const SPHParticles<d>& particles, const Array<Vector<real, d>>& arr, const KernelSPH& kernel)
	{
		return particles.Neighbor_Sum_Wrapper(Viscosity_Laplacian_Term_Func<d>(particles, arr, kernel));
	}
	template SingleIFuncV<2> Viscosity_Laplacian_Func<2>(const SPHParticles<2>&, const Array<Vector<real, 2>>&, const KernelSPH&);
	template SingleIFuncV<3> Viscosity_Laplacian_Func<3>(const SPHParticles<3>&, const Array<Vector<real, 3>>&, const KernelSPH&);


	template<int d>
	PairIFuncS Mass_Density_Term_Func(const Array<Vector<real, d> >& x, const Array<real>& m, const KernelSPH& kernel) {
		return [&](const int& i, const int& j) {
			return m[j] * kernel.Weight<d>(x[i] - x[j]);
		};
	}
	template<int d>
	SingleIFuncS Mass_Density_Func(const SPHParticles<d>& particles, const KernelSPH& kernel) {
		return particles.Neighbor_Sum_Wrapper(Mass_Density_Term_Func<d>(particles.XRef(), particles.MRef(), kernel));
	}
	template SingleIFuncS Mass_Density_Func<2>(const SPHParticles<2>&, const KernelSPH&);
	template SingleIFuncS Mass_Density_Func<3>(const SPHParticles<3>&, const KernelSPH&);

	template<int d>
	SingleIFuncS State_Equation(const SPHParticles<d>& particles, const real& rest_density, const real& kp, const int& kb)
	{
		return [&,rest_density,kp,kb](const int& i) {
			return kp * (MathFunc::Quick_Pow(particles.Rho(i) / rest_density, kb) - 1.0);
		};
	}
	template SingleIFuncS State_Equation<2>(const SPHParticles<2>&, const real&, const real&, const int&);
	template SingleIFuncS State_Equation<3>(const SPHParticles<3>&, const real&, const real&, const int&);
}

template<int d>
real PCISPHSolver<d>::Solve(SPHParticles<d>& particles, const real& dt)
{
	////Prepare data structures
	static Array<Vector<real, d> > predicted_v, predicted_x;
	static Array<real> pressure, predicted_rho;
	predicted_v = particles.VRef();
	predicted_x.resize(particles.Size());
	pressure.resize(particles.Size());
	AuxFunc::Fill(pressure, 0);//initial guess: pressure=0
	predicted_rho.resize(particles.Size());

	////save information in prior to boost calculation, because they will keep unchanged in this solver.
	//m_j * grad_w(i,j)
	static Array<Array<VectorD> > neighbor_mgrads; neighbor_mgrads.resize(particles.Size());
	particles.Exec_Each(
		[&](const int i) {
		const Array<int>& nbs = particles.nbs_searcher->Neighbors(i);
		size_t n = nbs.size();
		neighbor_mgrads[i].resize(n);
		for (size_t k = 0; k < n; k++) {
			int j = nbs[k];
			VectorD rij = particles.X(i) - particles.X(j);
			neighbor_mgrads[i][k] = particles.M(j) * grad_kernel.Grad<d>(rij);
		}
	}
	);

	//"state equation": p_i=beta*(rho_0-rho_i)
	real beta = alpha / (dt * dt);
	SingleIFuncS pred_rho_func = particles.Neighbor_Sum_Wrapper(SPHFunc::Mass_Density_Term_Func<d>(predicted_x, particles.MRef(), particles.kernel));
	SingleIFunc<void> update_p_func = [&](const int i) {
		pressure[i] += beta * (rho_0 - predicted_rho[i]);
		pressure[i] = std::max(pressure[i], 0.0);
	};
	//SingleIFuncV<d> gradp_func = SPHFunc::Symmetric_Grad_Func<d>(particles, pressure, grad_kernel);
	SingleIFuncV<d> gradp_func = [&](const int i) {
		const Array<int>& nbs = particles.nbs_searcher->Neighbors(i);
		size_t n = nbs.size();
		VectorD gradp = VectorD::Zero();
		const real& rho_i = particles.Rho(i);
		for (size_t k = 0; k < n; k++) {
			int j = nbs[k];
			const real& rho_j = particles.Rho(j);
			gradp += (pressure[i] / (rho_i * rho_i) + pressure[j] / (rho_j * rho_j)) * neighbor_mgrads[i][k];
		}
		return gradp * particles.Rho(i);
	};
	SingleIFuncV<d> pred_v_func = [&](const int i) {
		return particles.V(i) - dt / particles.Rho(i) * gradp_func(i);
	};
	
	int iter;
	real rho_err = 0;
	for (iter = 0; iter < max_iter; iter++) {
		//step 1: move particles and compute new density
		ArrayFunc::Copy(predicted_x, particles.XRef());
		intg_x_ptr->Integrate(predicted_x, predicted_v, dt);
		particles.Calc_Each(pred_rho_func, predicted_rho);
		rho_err = ArrayFunc::Upper_Rel_Error(rho_0, predicted_rho);
		if (rho_err <= max_rel_err) break;
		//step 2: use the predicted density to get a better guess for pressure
		particles.Exec_Each(update_p_func);
		//step 3: use gussed pressure to modify predicted velocity
		particles.Calc_Each(pred_v_func, predicted_v);
	}
	if (iter == max_iter) {
		std::cout << "[WARNING]PCISPH_Solver<d>::Solve exceeds " << max_iter << " iterations with error " << rho_err << "\n";
	}
	ArrayFunc::Copy(particles.VRef(), predicted_v);
	if (verbose) std::cout << "PCI carried " << iter << " steps with error " << rho_err << ", mean rho " << AuxFunc::Mean(particles.RhoRef()) << ", max rho " << AuxFunc::Max(particles.RhoRef()) << ", min rho" << AuxFunc::Min(particles.RhoRef()) << "\n";
	return rho_err;
}

template class PCISPHSolver<2>;
template class PCISPHSolver<3>;

template<int d>
real IISPHSolver<d>::Solve(SPHParticles<d>& particles, const real& dt)
{
	////Part I: figure out predicted density: rho_star without projection
	static Array<Vector<real, d> > predicted_v, predicted_x;
	predicted_v = particles.VRef(), predicted_x = particles.XRef();
	intg_x_ptr->Integrate(predicted_x, predicted_v, dt);
	static Array<real> rho_star; rho_star.resize(particles.Size());
	SingleIFuncS rho_star_func = particles.Neighbor_Sum_Wrapper(SPHFunc::Mass_Density_Term_Func<d>(predicted_x, particles.MRef(), particles.kernel));
	particles.Calc_Each(rho_star_func, rho_star);

	////Part II: figure out a_ii. They will be unchanged during the iterations.
	static Array<VectorD> grad_m; grad_m.resize(particles.Size());
	static Array<real> diag_terms; diag_terms.resize(particles.Size());
	//positions->grad_m.
	SingleIFuncV<d> grad_m_func = particles.template Neighbor_Sum_Wrapper<VectorD>(
		[&](const int i, const int j) {
		return particles.M(j) * particles.kernel.template Grad<d>(particles.X(i) - particles.X(j));
	}
	);
	//grad_m->diag_terms
	SingleIFuncS diag_term_func = (-dt * dt) * particles.template Neighbor_Sum_Wrapper<real>(
		[&](const int i, const int j) {
		VectorD Wij = particles.kernel.template Grad<d>(particles.X(i) - particles.X(j));
		return particles.M(j) / MathFunc::Power2(particles.Rho(i)) * Wij.dot(grad_m[i] + particles.M(i) * Wij);
	}
	);
	particles.Calc_Each(grad_m_func, grad_m);
	particles.Calc_Each(diag_term_func, diag_terms);

	////Part III: figure out LHS: (Ap)_i
	static Array<real> pressure; ArrayFunc::Resize_To<real>(pressure, particles.Size(), 0);//initial pressure guess: all 0
	static Array<VectorD> pressure_acc; pressure_acc.resize(particles.Size());
	//pressure->pressure acceleration, symmetric gradient form
	SingleIFuncV<d> p_acc_func = particles.template Neighbor_Sum_Wrapper<VectorD>(
		[&](const int i, const int j) {
		return -particles.M(j) *
			(pressure[i] / MathFunc::Power2(particles.Rho(i)) + pressure[j] / MathFunc::Power2(particles.Rho(j)))
			* particles.kernel.template Grad<d>(particles.X(i) - particles.X(j));
	}
	);
	//pressure acceleration->(Ap)_i
	//note that s_i=rho_0-rho_i
	SingleIFuncS LHS_func = (dt * dt) * particles.template Neighbor_Sum_Wrapper<real>(
		[&](const int i, const int j) {
		return particles.M(j) * (pressure_acc[i] - pressure_acc[j]).dot(particles.kernel.template Grad<d>(particles.X(i) - particles.X(j)));
	}
	);

	////Part IV: carry out the iteration. Note that diagonal elements and s_i=rho_0-rho_i will be unchanged
	static Array<real> predicted_rhos; predicted_rhos.resize(particles.Size());
	SingleIFunc<void> pressure_update_func = [&](const int i) {
		predicted_rhos[i] = rho_star[i] + LHS_func(i);
		real p = pressure[i] + omega / diag_terms[i] * (rho_0 - predicted_rhos[i]);
		pressure[i] = std::max(p, 0.0);
	};
	int iter; real rho_err = 0;
	for (iter = 0; iter < max_iter; iter++) {
		particles.Calc_Each(p_acc_func, pressure_acc);
		particles.Exec_Each(pressure_update_func);
		rho_err = ArrayFunc::Upper_Rel_Error(rho_0, predicted_rhos);
		if (rho_err < max_rel_err) break;
	}
	//update velocities with respect to pressure
	ArrayFunc::Array_Add(particles.VRef(), pressure_acc, dt);
	if (iter == max_iter) {
		std::cout << "[WARNING]IISPH_Solver<d>::Solve exceeds " << max_iter << " iterations with error " << rho_err << "\n";
	}
	if (verbose) std::cout << "#     IISPH carried " << iter << " steps with error " << rho_err << ", mean rho " << AuxFunc::Mean(particles.RhoRef()) << ", max rho " << AuxFunc::Max(particles.RhoRef()) << ", min rho" << AuxFunc::Min(particles.RhoRef()) << "\n";
	return real();
}
template class IISPHSolver<2>;
template class IISPHSolver<3>;
