//////////////////////////////////////////////////////////////////////////
// SPH Fluid
// Copyright (c) (2018-), Xiangxin Kong, Bo Zhu, Megndi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __FluidSPH_h__
#define __FluidSPH_h__
#include "Field.h"
#include "SPHFunc.h"
#include "SPHParticles.h"
#include "Params.h"
#include "Integrators.h"

class ParamsFluidSPH : public ParamsSPHParticles {
public:
	Declare_Param(real, dx);
	Declare_Param(Vector3, g);
	Declare_Param(real, band_dx_num);
	Register_Params(dx, g, band_dx_num);
	ParamsFluidSPH(const Vector3& _g, const real& _dx, const real& _kernel_dx_num = 3, const KernelType& mode = KernelType::SPIKY, const real& _band_dx_num = 1) :
		ParamsSPHParticles(_dx* _kernel_dx_num, mode)
	{
		dx = _dx;
		g = _g;
		band_dx_num = _band_dx_num;
		Register_Attributes();
	}
};

template<int d>
class FluidSPH {
	Typedef_VectorD(d);
public:
	SPHParticles<d> particles;
	real dx;
	real rest_number_density;
	VectorD gravity;//gravitational acceleration
	bool verbose = false;
	std::shared_ptr<IntegratorX<d>> intg_x_ptr;
	FluidSPH(const ParamsFluidSPH& params, std::shared_ptr<IntegratorX<d>> _intg)
		:particles(params),
		intg_x_ptr(_intg)
	{
		dx = params.dx;
		rest_number_density = SPHFunc::Rest_Number_Density<d>(dx, particles.kernel);
		gravity = AuxFunc::V<d>((Vector3)params.g);
	}
	void Update_Neighbors(void);
	virtual void Advance(const real dt, const real time = 0) = 0;
	real Max_Velocity(void);
};


class ParamsFluidSPHWC : public ParamsFluidSPH {
public:
	Declare_Param(real, kp);
	Declare_Param(real, kb);//state equation: p=kp*((rho/rho_0)^kb-1)
	Declare_Param(real, particle_mass);//all particles have the same mass
	Declare_Param(real, vis);//kinematic viscosity
	Register_Params(kp, kb, particle_mass, vis);
	ParamsFluidSPHWC(const ParamsFluidSPH& sph, const real& _kp, const real& _kb, const real &_mass, const real &_vis)
		:ParamsFluidSPH(sph)
	{
		kp = _kp;
		kb = _kb;
		particle_mass = _mass;
		vis = _vis;
		Register_Attributes();
	}
};

template<int d>
class FluidSPHWC : public FluidSPH<d> {
	Typedef_VectorD(d);
	using Base = FluidSPH<d>;
public:
	using Base::particles;
	real kp, kb;//state equation: p=kp*((rho/rho_0)^kb-1)
	real rest_mass_density;
	real vis_mu;//dynamic viscosity
	FluidSPHWC(const ParamsFluidSPHWC& params, const ImplicitShape<d>& boundary) :
		Base(params,
			std::make_shared<IntegratorXImplicit<d>>(boundary, 1.0, 0.0)
		)
	{
		kp = params.kp;
		kb = params.kb;
		rest_mass_density = Base::rest_number_density * params.particle_mass;
		vis_mu = params.vis * rest_mass_density;
	}
	virtual void Advance(const real dt, const real time = 0);
};

class ParamsFluidSPHPCI : public ParamsFluidSPH {
public:
	Declare_Param(real, particle_mass);//all particles have the same mass
	Declare_Param(real, vis);//kinematic viscosity
	Register_Params(particle_mass, vis);
	ParamsFluidSPHPCI(const ParamsFluidSPH& sph, const real _mass, const real _vis) :
		ParamsFluidSPH(sph)
	{
		particle_mass = _mass;
		vis = _vis;
		Register_Attributes();
	}
};

template<int d>
class FluidSPHPCI : public FluidSPH<d> {
	using Base = FluidSPH<d>;
public:
	using Base::particles;
	Typedef_VectorD(d);
	PCISPHSolver<d> pci_solver;
	real vis_mu;
	FluidSPHPCI(const ParamsFluidSPHPCI& params, const ImplicitShape<d> &boundary) :
		Base(params,
			std::make_shared<IntegratorXImplicit<d>>(boundary, 1.0, -0.1)
		),
		pci_solver(
			Base::intg_x_ptr,
			Base::rest_number_density* params.particle_mass,
			SPHFunc::Calculate_PCI_Coefficient<d>(Base::dx, Base::particles.kernel),
			KernelSPH(particles.kernel),
			10,//max iterations
			1e-2//max relative error
		)
	{
		vis_mu = params.vis * pci_solver.rho_0;
	}
	virtual void Advance(const real dt, const real time = 0);
};

class ParamsFluidIISPH : public ParamsFluidSPH {
public:
	Declare_Param(real, particle_mass);//all particles have the same mass
	Declare_Param(real, vis);//kinematic viscosity
	Register_Params(particle_mass, vis);
	ParamsFluidIISPH(const ParamsFluidSPH& sph, const real _mass, const real _vis) :
		ParamsFluidSPH(sph)
	{
		particle_mass = _mass;
		vis = _vis;
		Register_Attributes();
	}
};


template<int d>
class FluidIISPH : public FluidSPH<d> {
	using Base = FluidSPH<d>;
public:
	using Base::particles;
	Typedef_VectorD(d);
	IISPHSolver<d> solver;
	real vis_mu;
	FluidIISPH(const ParamsFluidIISPH& params, const ImplicitShape<d>& boundary) :
		Base(params,
			std::make_shared<IntegratorXImplicit<d>>(boundary, 1.0, -0.1)
		),
		solver(
			Base::intg_x_ptr,//integrator
			Base::rest_number_density* params.particle_mass,//rho_0
			0.05,//omega
			50,//max_iter
			1e-3//max_error
		) 
	{
		vis_mu = params.vis * solver.rho_0;
	}
	virtual void Advance(const real dt, const real time = 0);
};

//number density
template<int d,int bs=128> class FluidSPHND
{Typedef_VectorDii(d);
public:
	NDSPHParticles<d> particles;
	Array<ArrayF<int,bs*4> > nbs;
	SpatialHashing<d,bs> spatial_hashing;
	std::shared_ptr< KernelSPH > kernel;
	std::function<void(const real)> Collision;

	real length_scale;	////typical distance between two neighboring particles
	real den_0;			////rest density	
	real nden_0;		////rest number density
	real V_0;			////rest particle volume
	real mass_0;		////rest particle mass
	int avg_nb_num;		////averaged #particles
	real h;				////supporting radius			
	real kp;			////density-pressure coefficient
	real vis;			////viscosity
	real st;			////surface tension

	VectorD g=VectorD::Unit(1)*(real)-1.;
	bool use_body_force=false;
	bool use_fixed_dt=false;
	bool use_central_gravity=false;
	bool verbose=false;

public:
	void Initialize(const real _length_scale,const real _rho0=(real)1000.)
	{
		length_scale=_length_scale;
		den_0=_rho0;
		avg_nb_num=16*(int)pow(2,d-1);	////1d - 16; 2d - 32; 3d - 64

		if constexpr (d==1){
			V_0=length_scale;
			h=avg_nb_num*V_0;}
		else if constexpr (d==2){
			V_0=pi*pow((real).5*length_scale,2);
			h=sqrt((real)avg_nb_num*V_0/pi);}
		else if constexpr (d==3){
			V_0=(real)4/(real)3*pi*pow((real).5*length_scale,3);
			h=pow((real)3*(real)avg_nb_num*V_0/((real)4*pi),one_third);}

		mass_0=den_0*V_0;
		kp=(real)1e3;
		vis=den_0*(real)1e-2;
		st=(real)1;

		kernel=std::make_shared<KernelSPH>(h);
		spatial_hashing.Initialize(h);

		////update nden_0
		nbs.resize(particles.Size());
		Update_Neighbors();
		int pn=particles.Size();
		nden_0=(real)0;
		for(int i=0;i<pn;i++){
			real nden=(real)0;
			int nb_n=nbs[i][0];
			for (int k = 1; k <= nb_n; k++) {
				int j = nbs[i][k];
				real dis_ij = (particles.X(i) - particles.X(j)).norm();
				nden += kernel->Weight(d, dis_ij, KernelType::POLY6);
			}
			particles.M(i)=den_0/nden;
			nden_0+=nden;}
		nden_0/=(real)pn;

		if(verbose)Print_Params();
	}

	void Update_Neighbors()
	{
		spatial_hashing.Update_Voxels(particles.XRef());
		int pn=particles.Size();
		for(int i=0;i<pn;i++){
			spatial_hashing.Find_Nbs(particles.X(i),particles.XRef(),h,nbs[i]);}	
	}

	virtual void Advance(const real dt,const real time=0)
	{	
		////update nbs
		Update_Neighbors();
		const int pn=particles.Size();

		////update number density and pressure
		for(int i=0;i<pn;i++){
			real nden=(real)0;
			int nb_n=nbs[i][0];
			for (int k = 1; k <= nb_n; k++) {
				int j = nbs[i][k];
				real dis_ij = (particles.X(i) - particles.X(j)).norm();
				nden += kernel->Weight(d, dis_ij, KernelType::QUINTIC);
			}
			particles.ND(i)=nden;
			particles.Vol(i)=(real)1/nden;
			particles.P(i)=kp*(nden/nden_0-(real)1.);}

		////update forces
		for(int i=0;i<pn;i++){
			real one_over_m=(real)1/particles.M(i);
			VectorD f=VectorD::Zero();
			int nb_n=nbs[i][0];
			for(int k=1;k<=nb_n;k++){
				int j=nbs[i][k];
				VectorD r_ij=particles.X(i)-particles.X(j);
				real r2=r_ij.squaredNorm();
				real one_over_r2=(r2==(real)0?(real)0:(real)1/r2);
				VectorD v_ij=particles.V(i)-particles.V(j);
				VectorD f_p = -(particles.P(i) * pow(particles.Vol(i), 2) + particles.P(j) * pow(particles.Vol(j), 2)) * kernel->Grad<d>(r_ij, KernelType::SPIKY);
				VectorD f_v = vis * particles.Vol(i) * particles.Vol(j) * v_ij * one_over_r2 * r_ij.dot(kernel->Grad<d>(r_ij, KernelType::SPIKY));
				f+=(f_p+f_v);}
			VectorD a_g=use_central_gravity?-particles.X(i).normalized():g;
			particles.F(i)=one_over_m*f+a_g;}

		////time integration
		for(int i=0;i<pn;i++){
			particles.V(i)+=particles.F(i)*dt;
			particles.X(i)+=particles.V(i)*dt;}

		////collision
		if(Collision!=nullptr)Collision(dt);

		if(verbose)Print_Statistics();
	}

	real CFL() const
	{
		if(use_fixed_dt) return (real)1;
		else{
			real max_abs=(real)1e-5;
			for(int i=0;i<particles.Size();i++){
				real v_mag=particles.V(i).norm();
				if(v_mag>max_abs)max_abs=v_mag;}
			return length_scale/max_abs;}
	}

	void Print_Params()
	{
		std::cout<<"length_scale: "<<length_scale<<std::endl;
		std::cout<<"den_0: "<<den_0<<std::endl;
		std::cout<<"nden_0: "<<nden_0<<std::endl;
		std::cout<<"V_0: "<<V_0<<std::endl;
		std::cout<<"mass_0: "<<mass_0<<std::endl;
		std::cout<<"avg_nb_num: "<<avg_nb_num<<std::endl;
		std::cout<<"h: "<<h<<std::endl;
		std::cout<<"kp: "<<kp<<std::endl;
		std::cout<<"vis: "<<vis<<std::endl;
		std::cout<<"st: "<<st<<std::endl;
	}

	void Print_Statistics()
	{
		int pn=particles.Size();

		int avg_nb_num=0;
		int max_nb_num=-1;
		int min_nb_num=std::numeric_limits<int>::max();
		real avg_nden=(real)0;
		real max_nden=(real)-1;
		real min_nden=std::numeric_limits<real>::max();
		for(int i=0;i<pn;i++){
			avg_nb_num+=nbs[i][0];
			if(nbs[i][0]>max_nb_num)max_nb_num=nbs[i][0];
			if(nbs[i][0]<min_nb_num)min_nb_num=nbs[i][0];
			avg_nden+=particles.ND(i);
			if(particles.ND(i)>max_nden)max_nden=particles.ND(i);
			if(particles.ND(i)<min_nden)min_nden=particles.ND(i);}

		std::cout<<"nbs_num avg: "<<avg_nb_num/pn<<", max: "<<max_nb_num<<", min: "<<min_nb_num<<std::endl;
		std::cout<<"nden init: "<<nden_0<<", avg: "<<avg_nden/(real)pn<<", max: "<<max_nden<<", min: "<<min_nden<<std::endl;
	}

	//virtual void Update_Color_Field_And_Surface_Normal()
	//{
	//	for (int i=0; i<particles.Size(); i++)
	//	{
	//		real i_color=(real)0;
	//		VectorD i_sn=VectorD::Zero();
	//		for (int k=1; k<=nbs[i][0]; k++)
	//		{
	//			int j=nbs[i][k];
	//			VectorD r_ij=particles.X(i)-particles.X(j);
	//			real dist_ij=r_ij.norm();
	//			i_sn += particles.M(j)/particles.D(j)*kernel->Gradient_poly6(r_ij);
	//		}
	//		particles.C(i)=i_color;
	//		particles.SN(i)=i_sn;
	//	}
	//}

	//virtual void Update_Curvature_And_Surface_Force()
	//{
	//	for (int i=0; i<particles.Size(); i++){
	//		real i_curvature=(real)0;
	//		real surface_norm=particles.SN(i).norm();
	//		VectorD i_force=VectorD::Zero();
	//		for (int k=1; k<=nbs[i][0]; k++){
	//			int j=nbs[i][k];
	//			VectorD r_ij=particles.X(i)-particles.X(j);
	//			real dist_ij=r_ij.norm();				
	//			if(surface_norm != (real)0)
	//				i_curvature -= particles.M(j)/particles.D(j)*kernel->Wpoly6(dist_ij)/surface_norm;
	//		}
	//		particles.Curvature(i)=i_curvature;

	//		
	//		if (surface_norm > surface_threshold){
	//			real coef=surface_tension*i_curvature;
	//			i_force=coef *particles.SN(i);
	//			particles.F(i) -= i_force;
	//		}				
	//	}
	//}


	//void Enforce_Boundary_Conditions()
	//{
	//	for (int i=0; i<particles.Size(); i++)
	//	{
	//		for (int j=0; j<env_objects.size(); j++)
	//		{
	//			real phi=env_objects[j]->Phi(particles.X(i));
	//			if (phi<particles.R(i))
	//			{
	//				VectorD normal=env_objects[j]->Normal(particles.X(i));
	//				particles.F(i) += normal*kd*(particles.R(i)-phi)*particles.D(i);
	//			}
	//		}
	//	}
	//}
};
#endif

