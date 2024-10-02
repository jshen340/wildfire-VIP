//#####################################################################
// SPH Fluid driver
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//#####################################################################
#include "FluidSPHDriver.h"
#include "PointSetFunc.h"

template<int d>
void FluidSPHDriver<d>::Case_1()
{
	VectorDi cell_counts = VectorDi::Ones() * 16;
	real dx = (real)1 / (real)cell_counts[0];
	VectorD domain_min = -(real).5 * dx * cell_counts.template cast<real>();
	Grid<d> grid(cell_counts, dx, domain_min);
	int pn = grid.cell_counts.prod();
	fluid.particles.Resize(pn);
	iterate_cell(iter, grid) {
		const VectorDi& cell = iter.Coord();
		int c = grid.Cell_Index(cell);
		VectorD pos = grid.Center(cell);
		Initialize_Particle(c, pos);
	}
	fluid.Initialize(dx);
	fluid.use_central_gravity = true;
}

template<int d>
void FluidSPHDriver<d>::Case_2()
{
	//cfl=(real).1;
	VectorDi cell_counts = VectorDi::Ones() * scale; cell_counts[0] /= 2;
	real dx = (real)1 / (real)scale; VectorD domain_min = dx * VectorD::Ones();
	Grid<d> grid(cell_counts, dx, domain_min);
	int pn = grid.cell_counts.prod();
	fluid.particles.Resize(pn);
	iterate_cell(iter, grid) {
		const VectorDi& cell = iter.Coord();
		int c = grid.Cell_Index(cell);
		VectorD pos = grid.Center(cell);
		Initialize_Particle(c, pos);
	}
	fluid.Initialize(dx);
	fluid.Collision = std::bind(&FluidSPHDriver::Collision, this, std::placeholders::_1);
}

template<int d>
void FluidSPHDriver<d>::Case_3(void)
{
	VectorD center = VectorD::Zero();
	real dx = 1;
	real R = 3 * dx;
	fluid.particles.Resize(1);
	Initialize_Particle(0, center);
	if constexpr (d == 2) {
		PointSetFunc::Extend_Hexagon_Packing(fluid.particles, R, dx);
		for (int i = 0; i < 100; i++) {
			real alpha = i / 100.0 * 2 * pi;
			Vector2 rv = Vector2(cos(alpha), sin(alpha)) * R*1.01;
			int idx = fluid.particles.Add_Element();
			fluid.particles.Copy_Element_From(idx, fluid.particles, 0);
			fluid.particles.X(idx) = center + rv;
		}
	}
	else if constexpr (d == 3) {
		MacGrid<3> mac_grid = PointSetFunc::Extend_FCC_Packing(fluid.particles, R, dx);
		mac_grid.grid.Write_To_File_3d(output_dir + "/0/grid");
	}
	
	//real alpha = SPHFunc::Calculate_PCI_Coefficient<d>(dx, R, KernelType::SPIKY);

	//NeighborKDTree<d> nbs_searcher;
	//nbs_searcher.Build_Data(fluid.particles.XRef());
	//Array<int> nbs;
	//nbs_searcher.Find_K_Nearest_Nb(fluid.particles.X(0), 100, nbs);
	//for (int i = 0; i < nbs.size(); i++) {
	//	int j = nbs[i];
	//	if (j != 0) {
	//		std::cout << fluid.particles.X(j).norm() << " ";
	//	}
	//}

	fluid.Initialize(dx);
}


template class FluidSPHDriver<2>;
template class FluidSPHDriver<3>;