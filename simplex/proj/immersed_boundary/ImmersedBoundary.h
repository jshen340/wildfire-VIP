//////////////////////////////////////////////////////////////////////////
// Copyright (c) (2021-), Fan Feng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __ImmersedBoundary_h__
#define __ImmersedBoundary_h__
#include "SoftBodyMassSpring.h"
#include "FluidEuler.h"
#include "Interpolation.h"
#include "NeighborSearcher.h"
#include "Particles.h"

template<int d> class ImmersedBoundary
{Typedef_VectorDii(d);
public:
	SoftBodyMassSpring<d> soft_body;
	FluidEuler<d> fluid;

	void Initialize(const VectorDi& cell_counts, const real dx, const VectorD& domain_min = VectorD::Zero())
	{
		fluid.Initialize(cell_counts,dx,domain_min);
		std::cout << "cell_counts: "<< cell_counts.transpose() << std::endl;
		std::cout << "dx: "<< dx << std::endl;
	}

	void Advance(const real dt)
	{
		Advect_Soft_Body(dt);
		Add_Force_To_Fluid_f(dt);
		fluid.Enforce_Incompressibility();
		Soft_Body_Collision();
		fluid.Advection(dt);
	}

	void Advect_Soft_Body(const real dt)
	{
		Interpolation<d> intp(fluid.mac_grid);
		#pragma omp parallel for
		for (int i = 0; i < soft_body.particles.Size(); i++) {
			VectorD& pos = soft_body.particles.X(i);
			VectorD& v=intp.Interpolate_Face_Vectors_Quadratic(fluid.velocity, pos);
			pos += v * dt;
		}
	}

	void Add_Force_To_Fluid_f(real dt)
	{
		FaceField<real, d> v_w(fluid.mac_grid.grid.cell_counts, (real)0);
		FaceField<real, d> delta_v(fluid.mac_grid.grid.cell_counts, (real)0);

		Particles<d>& particles = soft_body.particles;
		int pn = particles.Size();
		Interpolation<d> intp(fluid.mac_grid.grid);

		soft_body.Update_Forces(particles.XRef(), particles.VRef(), particles.FRef());

		for (int i = 0; i < pn; i++) {
			const VectorD& pos = particles.X(i);
			VectorD vel_d = particles.F(i) * dt;
			intp.Interpolate_Point_To_Faces_Quadratic(pos, vel_d, delta_v, v_w);
		}

		////Add force from soft_body to fluid
		for (int axis = 0; axis < d; axis++) {
			int face_num = fluid.mac_grid.Number_Of_Faces(axis);
			#pragma omp parallel for
			for (int i = 0; i < face_num; i++) {
				VectorDi face = fluid.mac_grid.Face_Coord(axis, i);
				if (v_w(axis, face) != (real)0) {
					fluid.velocity(axis, face) += delta_v(axis, face) / v_w(axis, face);
				}
			}
		}
	}

	void Soft_Body_Collision() {
		int pn= soft_body.particles.Size();
		for (int i = 0; i < pn; i++) {
			VectorD& pos = soft_body.particles.X(i);
			for (int j = 0; j < d; j++) {
				if (pos[j] < fluid.mac_grid.grid.domain_min[j])pos[j] = fluid.mac_grid.grid.domain_min[j];
				else if (pos[j] > fluid.mac_grid.grid.domain_max[j])pos[j] = fluid.mac_grid.grid.domain_max[j];
			}
		}
	}
};
#endif
