//////////////////////////////////////////////////////////////////////////
// Project a vector field with free boundary to divergence free on a MAC grid
// Copyright (c) (2018-),Fan Feng, Shuqi Yang, Bo Zhu
// This file is part of SimpleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifdef USE_CPX


#include "ProjectionIrregularCPX.h"
#include "Timer.h"

//////////////////////////////////////////////////////////////////////////
////projection functions
template<class T, int d> void ProjectionIrregularCPX<T, d>::Apply_Jump_Condition_To_b(const T dt)
{
	T one_over_dx = 1 / mac_grid->grid.dx;

	int cell_num=mac_grid->grid.cell_counts.prod();
	#pragma omp parallel for
	for (int i = 0; i < cell_num; i++) {
		VectorDi cell = mac_grid->grid.Cell_Coord(i);
		if (!Is_Fluid_Cell(cell))continue;

		for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
			int axis;VectorDi face;MacGrid<d>::Cell_Incident_Face(cell,i,axis,face);
			real coef=!(*bc).Is_Psi_N(axis,face)?(real)1:(real)0;

			T theta; bool is_intf; T intf_coef = Intf_Coef(cell, i, theta, is_intf);
			if (is_intf) {
				VectorDi nb_cell = Grid<d>::Nb_C(cell, i);
				Vector<T, d> intf_pos = (1.0 - theta) * mac_grid->grid.Center(cell) + theta * mac_grid->grid.Center(nb_cell);
				T p_jump = Pressure_Jump(dt, intf_pos);
				cell_b(cell) += coef * intf_coef * p_jump * one_over_dx;
			}
		}
	}
}

////need to specify target_vol and current_vol 
template<class T, int d> void ProjectionIrregularCPX<T, d>::Apply_Vol_Control_To_b()
{
	if(!use_vol_control)return;
	Assert(target_vol > 0, "ProjectionIrregularCPX: target_vol not set");

	if (calc_current_vol)current_vol = levelset->Total_Volume();
	T vol_correction=(target_vol-current_vol)/current_vol;
	if (verbose)Info("vol correction: {}", vol_correction);

	int cell_num=mac_grid->grid.cell_counts.prod();
	#pragma omp parallel for
	for (int i = 0; i < cell_num; i++) {
		VectorDi cell = mac_grid->grid.Cell_Coord(i);
		if (!Is_Fluid_Cell(cell))continue;
		T cell_div = vol_correction;
		cell_b(cell) += vol_control_ks * mac_grid->grid.dx * cell_div;
	}
}

template<class T, int d> void ProjectionIrregularCPX<T, d>::Build(const T dt, const FaceField<T, d>& velocity)
{
	Timer timer;					timer.Reset();
	Base::Build_A(dt);				if(verbose)timer.Elapse_And_Output_And_Reset("Build A");
	Base::Build_b(dt,velocity);		if(verbose)timer.Elapse_And_Output_And_Reset("Build b");
	Apply_Jump_Condition_To_b(dt);	if(verbose)timer.Elapse_And_Output_And_Reset("Apply jump condition to b");
	Apply_Vol_Control_To_b();		if(verbose)timer.Elapse_And_Output_And_Reset("Apply volume control to b");
}

template<class T, int d> T ProjectionIrregularCPX<T, d>::Intf_Coef(const VectorDi& cell,const int i,T& theta,bool& is_intf) const
{
	VectorDi nb_cell=Grid<d>::Nb_C(cell,i);
	if(mac_grid->grid.Valid_Cell(nb_cell)){
		T phi0=(*levelset).phi(cell);T phi1=(*levelset).phi(nb_cell);
		return Intf_Coef(phi0,phi1,theta,is_intf);}
	return 1;
}

template<class T, int d> T ProjectionIrregularCPX<T, d>::Intf_Coef(const T& phi0, const T& phi1, T& theta, bool& is_intf) const
{
	T max_scale = 1e3; is_intf = false;
	if (LevelSet<d>::Interface(phi0, phi1)) {
		is_intf = true;
		theta = LevelSet<d>::Theta(phi0, phi1);
		T scale = std::min(1 / theta, max_scale); return scale;
	}
	return 1;
}

//////////////////////////////////////////////////////////////////////////
////Physical interface functions
//see: https://wmdcstdio.com/2021/07/11/projection-matrix-terms/
template<class T, int d> T ProjectionIrregularCPX<T, d>::Off_Diag_Term(const T dt, const VectorDi& fluid_cell, const int& nbidx) const
{
	VectorDi nb_cell = Grid<d>::Nb_C(fluid_cell, nbidx);
	int axis; VectorDi face; MacGrid<d>::Cell_Incident_Face(fluid_cell, nbidx, axis, face);
	if (bc->Is_Psi_N(axis, face))return 0;
	if (Is_Fluid_Cell(nb_cell)) return -1;
	return 0;
}

template<class T, int d> T ProjectionIrregularCPX<T, d>::Diag_Face_Term(const T dt, const int& axis, const VectorDi& face) const
{
	VectorDi cell[2]; for (int i = 0; i < 2; i++)cell[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
	if (!Is_Fluid_Cell(cell[0])) std::swap(cell[0], cell[1]);
	if (!Is_Fluid_Cell(cell[0])) return 0;
	T fluid_phi = levelset->phi(cell[0]);
	if (bc->Is_Psi_N(axis, face)) return 0;
	if (mac_grid->grid.Valid_Cell(cell[1])) {
		T nb_phi = levelset->phi(cell[1]);
		T theta; bool is_intf; T intf_coef = Intf_Coef(fluid_phi, nb_phi, theta, is_intf);
		if (is_intf) return intf_coef;
	}
	return 1;
}

template<class T, int d> T ProjectionIrregularCPX<T, d>::Velocity_Correction(const T dt, const int& axis, const VectorDi& face) const
{
	if ((*bc).Is_Psi_N(axis, face)) { return 0; }

	VectorDi cell[2]; for (int i = 0; i < 2; i++)cell[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
	if (!Is_Fluid_Cell(cell[0]) && !Is_Fluid_Cell(cell[1]))return 0.;

	T cell_p[2]; for (int i = 0; i < 2; i++)cell_p[i] = Is_Fluid_Cell(cell[i]) ? cpx_poisson.x(cell[i]) : 0;		////air cell has zero pressure

	if (mac_grid->grid.Valid_Cell(cell[0]) && mac_grid->grid.Valid_Cell(cell[1])) {
		T one_over_dx = 1.0 / mac_grid->grid.dx;
		T phi[2]; for (int i = 0; i < 2; i++)phi[i] = (*levelset).phi(cell[i]);
		int f_idx = (Is_Fluid_Cell(cell[0]) ? 0 : 1); int a_idx = (f_idx == 0 ? 1 : 0);
		T theta; bool is_intf; T intf_coef = Intf_Coef(phi[f_idx], phi[a_idx], theta, is_intf);
		if (is_intf) {		////set the boundary cell pressure to be the pressure jump value multiplying intf coef
			cell_p[f_idx] *= intf_coef;
			if (apply_jump_condition) {
				VectorD intf_pos = (1 - theta) * mac_grid->grid.Center(cell[f_idx]) + theta * mac_grid->grid.Center(cell[a_idx]);
				T p_jump = Pressure_Jump(dt, intf_pos);
				cell_p[a_idx] = p_jump * one_over_dx * intf_coef;
			}
		}
	}

	return -(cell_p[1] - cell_p[0]);
}

template class ProjectionIrregularCPX<Scalar, 2>;
template class ProjectionIrregularCPX<Scalar, 3>;

#endif