//#####################################################################
// Universal SPH Driver
// Copyright (c) (2018-), Xiangxin Kong, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//#####################################################################
#include "SPHDriver.h"
#include "PointSetFunc.h"


template<int d>
void SPHDriver<d>::Add_SPH_Particle(SPHParticles<d>& particles, const VectorD& pos, const real& m)
{
	int idx = particles.Add_Element();
	particles.X(idx) = pos;
	particles.M(idx) = m;
	particles.V(idx) = VectorD::Zero();
	particles.Rho(idx) = 0;
}

template<int d>
void SPHDriver<d>::Case_1(void)
{
	VectorD center = VectorD::Zero();
	real bowl_radius = 2.0;
	real drop_radius = 1.0;
	//boundary
	auto boundary = ImplicitShape<d>::Bounding_Sphere(center, bowl_radius);
	//initialize solver
	real dx = 1.0 / scale;
	ParamsFluidSPH params_sph(Vector3(0, -1, 0), dx, 4, KernelType::SPIKY, 1.0);
	real mass = 1.0, vis = 1e-2;
	if (algo_type == 0) {
		real kp = 1e3, kb = 3;
		ParamsFluidSPHWC params_wc(params_sph, kp, kb, mass, vis);
		fluid = std::make_shared<FluidSPHWC<d>>(params_wc, boundary);
	}
	else if (algo_type == 1) {
		ParamsFluidSPHPCI params_pci(params_sph, mass, vis);
		fluid = std::make_shared<FluidSPHPCI<d>>(params_pci, boundary);
	}
	else if (algo_type == 2) {
		ParamsFluidIISPH params_ii(params_sph, mass, vis);
		fluid = std::make_shared<FluidIISPH<d>>(params_ii, boundary);
	}
	//add points
	Add_SPH_Particle(fluid->particles, center, mass);
	PointSetFunc::Extend_Inside<d>(fluid->particles, dx, drop_radius, nullptr);
}

template<int d>
void SPHDriver<d>::Case_2(void)
{
	VectorD center = VectorD::Zero();
	real bound_length = 2.0 * PhysicalUnits::m;
	real drop_length = 1.0 * PhysicalUnits::m;
	real drop_envelop_radius = drop_length / 2 * sqrt(3);
	real dx = drop_length / scale;
	//boundary and drop geometry
	ImplicitShape<d> boundary = ImplicitShape<d>::Bounding_Box(center, bound_length * VectorD::Ones());
	auto drop_ptr = std::make_shared<Box<d>>(center, drop_length);
	Points<d> points_positions;
	points_positions.Add_Position_Only(center);
	PointSetFunc::Extend_Inside(points_positions, dx, drop_envelop_radius, drop_ptr.get());
	//calculate particle mass and set all particles
	real water_volume = MathFunc::Quick_Pow(drop_length, d);
	real water_mass = water_volume * PhysicalConstants::rho_water;
	real mass = water_mass / points_positions.Size();//particle mass
	//initialize solver
	ParamsFluidSPH params_sph(Vector3(0, -1, 0) * PhysicalConstants::g, dx, 2, KernelType::CUBIC, 1.0);
	if (algo_type == 0) {
		real kp = 1e4, kb = 1, vis = 1e-2;
		ParamsFluidSPHWC params_wc(params_sph, kp, kb, mass, vis);
		fluid = std::make_shared<FluidSPHWC<d>>(params_wc, boundary);
	}
	else if (algo_type == 1) {
		ParamsFluidSPHPCI params_pci(params_sph, mass, PhysicalConstants::nu_water);
		fluid = std::make_shared<FluidSPHPCI<d>>(params_pci, boundary);
	}
	else if (algo_type == 2) {
		ParamsFluidIISPH params_ii(params_sph, mass, PhysicalConstants::nu_water);
		fluid = std::make_shared<FluidIISPH<d>>(params_ii, boundary);
	}
	Add_SPH_Particle(fluid->particles, center, mass);
	PointSetFunc::Extend_Identical_Particles<d>(fluid->particles, points_positions.XRef());
}

template<int d>
void SPHDriver<d>::Case_3(void)
{
	VectorD center = VectorD::Zero();
	real bound_length = 2.0 * PhysicalUnits::m;
	real dx = bound_length / scale;
	//boundary and drop geometry
	ImplicitShape<d> boundary = ImplicitShape<d>::Bounding_Box(center, bound_length * VectorD::Ones());
	VectorD drop_length = VectorD::Ones() * bound_length; drop_length[0] /= 2;
	VectorD drop_center = center - VectorD::Unit(0) * 0.25 * bound_length;
	VectorD drop_min = drop_center - drop_length / 2, drop_max = drop_center + drop_length / 2;
	auto drop_ptr = std::make_shared<Box<d>>(drop_min, drop_max);
	Points<d> points_positions;
	points_positions.Add_Position_Only(drop_center);
	PointSetFunc::Extend_Inside(points_positions, dx, drop_length.norm() / 2, drop_ptr.get());
	//calculate particle mass and set all particles
	real water_volume = drop_length.prod();
	real water_mass = water_volume * PhysicalConstants::rho_water;
	real mass = water_mass / points_positions.Size();//particle mass
	//initialize solver
	ParamsFluidSPH params_sph(Vector3(0, -1, 0) * PhysicalConstants::g, dx, 3, KernelType::CUBIC, 1.0);
	if (algo_type == 0) {
		real kp = 1e3, kb = 1, vis = 1e-2;
		ParamsFluidSPHWC params_wc(params_sph, kp, kb, mass, vis);
		fluid = std::make_shared<FluidSPHWC<d>>(params_wc, boundary);
	}
	else if (algo_type == 1) {
		ParamsFluidSPHPCI params_pci(params_sph, mass, PhysicalConstants::nu_water);
		fluid = std::make_shared<FluidSPHPCI<d>>(params_pci, boundary);
	}
	else if (algo_type == 2) {
		ParamsFluidIISPH params_ii(params_sph, mass, PhysicalConstants::nu_water);
		fluid = std::make_shared<FluidIISPH<d>>(params_ii, boundary);
	}
	Add_SPH_Particle(fluid->particles, center, mass);
	PointSetFunc::Extend_Identical_Particles<d>(fluid->particles, points_positions.XRef());
}

template<int d>
void SPHDriver<d>::Case_4(void)
{
	VectorD center = VectorD::Zero();
	real bound_length = 2.0 * PhysicalUnits::m;
	real dx = bound_length / scale;
	//boundary and drop geometry
	ImplicitShape<d> boundary = ImplicitShape<d>::Bounding_Box(center, bound_length * VectorD::Ones());
	VectorD drop_length = VectorD::Ones() * bound_length; drop_length[1] /= 2;
	VectorD drop_center = center - VectorD::Unit(1) * 0.25 * bound_length;
	VectorD drop_min = drop_center - drop_length / 2, drop_max = drop_center + drop_length / 2;
	auto drop_ptr = std::make_shared<Box<d>>(drop_min, drop_max);
	Points<d> points_positions;
	points_positions.Add_Position_Only(drop_center);
	PointSetFunc::Extend_Inside(points_positions, dx, drop_length.norm() / 2, drop_ptr.get());
	//calculate particle mass and set all particles
	real water_volume = drop_length.prod();
	real water_mass = water_volume * PhysicalConstants::rho_water;
	real mass = water_mass / points_positions.Size();//particle mass
	//initialize solver
	ParamsFluidSPH params_sph(Vector3(0, -1, 0) * PhysicalConstants::g, dx, 2, KernelType::CUBIC, 1.0);
	if (algo_type == 0) {
		real kp = 1e3, kb = 1, vis = 1e-2;
		ParamsFluidSPHWC params_wc(params_sph, kp, kb, mass, vis);
		fluid = std::make_shared<FluidSPHWC<d>>(params_wc, boundary);
	}
	else if (algo_type == 1) {
		ParamsFluidSPHPCI params_pci(params_sph, mass, PhysicalConstants::nu_water);
		fluid = std::make_shared<FluidSPHPCI<d>>(params_pci, boundary);
	}
	else if (algo_type == 2) {
		ParamsFluidIISPH params_ii(params_sph, mass, PhysicalConstants::nu_water);
		fluid = std::make_shared<FluidIISPH<d>>(params_ii, boundary);
	}
	Add_SPH_Particle(fluid->particles, center, mass);
	PointSetFunc::Extend_Identical_Particles<d>(fluid->particles, points_positions.XRef());
}

template class SPHDriver<2>;
template class SPHDriver<3>;