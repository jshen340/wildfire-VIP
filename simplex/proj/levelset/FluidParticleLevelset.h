#pragma once
#include <numeric>
#include "MacGrid.h"
#include "Advection.h"
#include "LevelSet.h"
#include "Timer.h"
#include "FluidEulerFreeSurface.h"
#include "ParticleLevelSet.h"

template<int d> class FluidParticleLevelset : public FluidEulerFreeSurface<d> {
    Typedef_VectorDii(d);
    using Base = FluidEulerFreeSurface<d>;
public:
    using Base::mac_grid;using Base::velocity;using Base::alpha;using Base::type;using Base::bc;using Base::g;using Base::use_body_force;
    using Base::levelset; using Base::narrow_band_width;
	using Base::Is_Fluid_Cell; using Base::Is_Fluid_Face;

    ParticleLevelSet<d> pls;
	bool use_particle=true;
	int adv_cnt=0;
	int reseed_cnt=20;

    virtual void Initialize_Interface()
	{
		Base::Initialize_Interface();
		pls.Initialize(mac_grid.grid);////Set band width here
	}

    virtual void Advection(const real dt)
	{
		////Advect velocity
		MacGrid<d> ghost_grid=mac_grid;FaceField<real,d> ghost_velocity=velocity;
		Advection::Semi_Lagrangian(dt,ghost_grid,ghost_velocity,mac_grid,velocity);
		//// ow
		if(use_particle){
			////Advect interface
			Advection_Single_Particle_Levelset(ghost_velocity, pls, dt);
			iterate_cell(iter, mac_grid.grid) { const auto& cell = iter.Coord();
					levelset.phi(cell)=pls.phi(cell);}
			adv_cnt++;
			if(adv_cnt%reseed_cnt==0) pls.Reseed();
			pls.Update_Particle_Radii();
		}else{
            ////Advect interface
            Field<real,d> ghost_phi=levelset.phi;
            // Advection::Semi_Lagrangian(dt,mac_grid,velocity,ghost_phi,mac_grid,levelset.phi);
			Advection::MacCormack(dt,mac_grid,velocity,levelset.phi);
            levelset.Fast_Marching();
		}


		Update_Cell_Types();
	}

    void Setup_Particle_Level_Set() {
        iterate_cell(iter,pls.grid){const VectorDi& cell=iter.Coord();
            pls.phi(cell)=levelset.phi(cell);}
        pls.Initialize_Particles();
        pls.Fast_Marching();
	}


    void Advection_Single_Particle_Levelset(const FaceField<real,d>& velocity, ParticleLevelSet<d>& pls, const real dt){
		Field<real,d> ghost_phi=pls.phi;
		// Advection::Semi_Lagrangian(dt,ghost_grid,ghost_velocity,ghost_phi,mac_grid,levelset.phi);
		Advection::MacCormack(dt,mac_grid,velocity,pls.phi);
		// pls.Fast_Marching();
		Interpolation<d> intp(pls.grid);

		////advect particles
		const bool pls_adv_rk2 = false;
		if(pls_adv_rk2){
			//// RK-2
			for(int i=0;i<pls.particles.Size();i++){
				VectorD v=intp.Interpolate_Face_Vectors(velocity, pls.particles.X(i));
				VectorD p2=pls.particles.X(i)+v*dt*(real).5;
				v=intp.Interpolate_Face_Vectors(velocity, p2);
				pls.particles.X(i)+=v*dt;}}
		else{
			//// RK-3
			const real one_ninth =1.0/9.0;
			for(int i=0;i<pls.particles.Size();i++){
				VectorD v1=intp.Interpolate_Face_Vectors(velocity, pls.particles.X(i));
				VectorD p2=pls.particles.X(i)+v1*dt*(real).5;
				VectorD v2=intp.Interpolate_Face_Vectors(velocity, p2);
				VectorD p3=pls.particles.X(i)+v2*dt*(real).75;
				VectorD v3=intp.Interpolate_Face_Vectors(velocity, p3);
				pls.particles.X(i)+=(v1*2+v2*3+v3*4)*dt*one_ninth;}}
		

		pls.Reseed_Escaped_Particles();
		pls.Correction();
		pls.Fast_Marching();////Reinitialization
		pls.Correction();
		pls.Update_Particle_Radii();
	}

};