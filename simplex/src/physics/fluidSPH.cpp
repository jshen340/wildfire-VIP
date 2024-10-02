#include "FluidSPH.h"
#include "ArrayFunc.h"

template<int d>
void FluidSPH<d>::Update_Neighbors(void)
{
	particles.Update();
}

template<int d>
real FluidSPH<d>::Max_Velocity(void)
{
	return AuxFunc::Linf_Norm(particles.VRef());
}

template class FluidSPH<2>;
template class FluidSPH<3>;

template<int d>
void FluidSPHWC<d>::Advance(const real dt, const real time)
{
	//initialize
	static Array<real> pressures;
	static Array<VectorD> forces;
	pressures.resize(particles.Size());
	forces.resize(particles.Size());
	AuxFunc::Fill(forces, VectorD::Zero());
	particles.Update();
	//pressure force
	SingleIFuncS rho_func = SPHFunc::Mass_Density_Func(particles, particles.kernel);
	particles.Calc_Each(rho_func, particles.RhoRef());
	SingleIFuncS p_func = SPHFunc::State_Equation(particles, rest_mass_density, kp, kb);
	particles.Calc_Each(p_func, pressures);
	SingleIFuncV<d> p_force_func = SPHFunc::Symmetric_Grad_Func(particles, pressures, particles.kernel) * (-1.0);
	particles.Add_Each(p_force_func, forces);
	//gravitational force
	SingleIFuncV<d> g_force_func = SPHFunc::Index_Function(particles.RhoRef()) * Base::gravity;
	particles.Add_Each(g_force_func, forces);
	//viscosity
	SingleIFuncV<d> vis_force_func = vis_mu * SPHFunc::Viscosity_Laplacian_Func<d>(particles, particles.VRef(), particles.kernel);
	particles.Add_Each(vis_force_func, forces);
	//time integration
	auto intg_v = [&](const int& i) {
		Vector<real, d> dv = dt * forces[i] / particles.Rho(i);
		particles.V(i) += dv;
	};
	particles.Exec_Each(intg_v);
	Base::intg_x_ptr->Integrate(particles.XRef(), particles.VRef(), dt);
}

template class FluidSPHWC<2>;
template class FluidSPHWC<3>;

template<int d>
void FluidSPHPCI<d>::Advance(const real dt, const real time)
{
	//initialize
	static Array<VectorD> forces;
	forces.resize(particles.Size());
	AuxFunc::Fill(forces, VectorD::Zero());
	particles.Update();
	//rho for dynamics, and external forces
	particles.template Calc_Each<real>(SPHFunc::Mass_Density_Func(particles, particles.kernel), particles.RhoRef());
	particles.template Add_Each<VectorD>(SPHFunc::Index_Function(particles.RhoRef()) * Base::gravity, forces);
	particles.template Add_Each<VectorD>(SPHFunc::Viscosity_Laplacian_Func<d>(particles, particles.VRef(), particles.kernel) * vis_mu, forces);
	//time integration
	auto intg_v = [&](const int& i) {
		Vector<real, d> dv = dt * forces[i] / particles.Rho(i);
		particles.V(i) += dv;
	};
	particles.Exec_Each(intg_v);
	pci_solver.Solve(particles, dt);
	Base::intg_x_ptr->Integrate(particles.XRef(), particles.VRef(), dt);
	ArrayFunc::Numerical_Check(particles.XRef(), "X", true);
}

template class FluidSPHPCI<2>;
template class FluidSPHPCI<3>;

template<int d>
void FluidIISPH<d>::Advance(const real dt, const real time)
{
	//initialize
	static Array<VectorD> forces;
	forces.resize(particles.Size());
	AuxFunc::Fill(forces, VectorD::Zero());
	particles.Update();
	//rho for dynamics, and external forces
	particles.template Calc_Each<real>(SPHFunc::Mass_Density_Func(particles, particles.kernel), particles.RhoRef());
	particles.template Add_Each<VectorD>(SPHFunc::Index_Function(particles.RhoRef()) * Base::gravity, forces);
	particles.template Add_Each<VectorD>(SPHFunc::Viscosity_Laplacian_Func<d>(particles, particles.VRef(), particles.kernel) * vis_mu, forces);
	//time integration
	auto intg_v = [&](const int& i) {
		Vector<real, d> dv = dt * forces[i] / particles.Rho(i);
		particles.V(i) += dv;
	};
	particles.Exec_Each(intg_v);
	solver.Solve(particles, dt);
	Base::intg_x_ptr->Integrate(particles.XRef(), particles.VRef(), dt);
	ArrayFunc::Numerical_Check(particles.XRef(), "X", true);
}


template class FluidIISPH<2>;
template class FluidIISPH<3>;
