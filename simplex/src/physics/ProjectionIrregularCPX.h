//////////////////////////////////////////////////////////////////////////
// Project a vector field with free boundary to divergence free on a MAC grid
// Copyright (c) (2018-), Fan Feng, Shuqi Yang, Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __ProjectionIrregularCPX_h__
#define __ProjectionIrregularCPX_h__

#ifdef USE_CPX

#include "ProjectionCPX.h"
#include "LevelSet.h"

template<class T, int d> class ProjectionIrregularCPX : public ProjectionCPX<T, d>
{Typedef_VectorDii(d);
public:
	using Base=ProjectionCPX<T, d>;
	using Base::mac_grid;
	using Base::bc;
	using Base::verbose;
	using Base::cell_b;
	using Base::cpx_poisson;

	const LevelSet<d>* levelset=nullptr;

	SurfaceTensionMode surface_tension_mode;

	//bool use_surface_tension=false;

	////jump condition (surface tension)
	//bool use_jump_condition=false;
	T sigma=(T)1e-5;					////surface tension coefficient for the default pressure jump

	////divergence control
	bool use_vol_control=false;
	bool calc_current_vol=true;						////calculate the current vol within Apply_Vol_Control_To_b or not
	T target_vol=(T)-1;						////ALWAYS NEEDS TO BE SET EXTERNALLY!
	T current_vol=(T)-1;						////needs to be set externally or calculated within Apply_Vol_Control_To_b when calc_current_vol=true			
	T vol_control_ks=(T)1e2;

	// default narrow band width
	int narrow_band_cell_num=5;

	//implicit gravity
	bool use_jump_gravity = false;
	VectorD gravity_coefficient = VectorD::Zero();//(rho_L-rho_A) * g

protected:
	bool apply_jump_condition = false;

	////callback functions
	//std::function<real(const VectorD&)> Jump_Condition = std::bind(&ProjectionIrregularCPX<T, d>::Pressure_Jump, this, std::placeholders::_1);

public:
	////constructors
	ProjectionIrregularCPX(json& j, const MacGrid<d>* _mac_grid, const LevelSet<d>* _levelset, const Field<ushort, d>* _type = nullptr, const BoundaryConditionMacGrid<d>* _bc = nullptr)
		:Base(j, _mac_grid, _type, _bc)
	{
		levelset = _levelset;
		surface_tension_mode = j.value<SurfaceTensionMode>("surface_tension_mode", SurfaceTensionMode::NONE);
		sigma = j.value<T>("sigma", 1e-5);
		use_vol_control = j.value<bool>("use_vol_control", false);
		calc_current_vol = j.value<bool>("calc_current_vol", true);
		target_vol = j.value<T>("target_vol", -1);
		current_vol = j.value<T>("current_vol", -1);
		vol_control_ks = j.value<T>("vol_control_ks", 1e2);
		narrow_band_cell_num = j.value<int>("narrow_band_cell_num", 5);
		use_jump_gravity = j.value<bool>("use_jump_gravity", false);
		gravity_coefficient = j.value<VectorD>("gravity_coefficient", VectorD::Zero());

		//calculate values
		apply_jump_condition = (surface_tension_mode == SurfaceTensionMode::EXPLICIT) || use_jump_gravity;
	}
	~ProjectionIrregularCPX() {

	}

	////projection functions
	virtual void Apply_Jump_Condition_To_b(const T dt);
	virtual void Apply_Vol_Control_To_b();
	virtual void Build(const T dt, const FaceField<T, d>& velocity);

	////interface coef returns clamped 1/theta
	T Intf_Coef(const T& phi0, const T& phi1, T& theta, bool& is_intf) const;				////levelset-related
	T Intf_Coef(const VectorDi& cell,/*nb idx*/const int i, T& theta, bool& is_intf) const;		////levelset-related

	//////////////////////////////////////////////////////////////////////////
	////default callback functions: all levelset-related

	inline virtual T Pressure_Jump(const T dt, const Vector<T, d>& pos) const {
		T p_jump = 0;
		if (surface_tension_mode == SurfaceTensionMode::EXPLICIT) {
			T curvature = (*levelset).Curvature(pos);
			p_jump += sigma * curvature;
		}
		if (use_jump_gravity) {
			p_jump -= gravity_coefficient.dot(pos);
		}
		return p_jump * dt;
	}

	inline bool Is_Levelset_Interface_Face_Index(const std::pair<int, int> face_idx) const
	{
		int axis=face_idx.first;
		VectorDi face=mac_grid->face_grids[axis].Node_Coord(face_idx.second);
		Vector<T, d> pos = mac_grid->Face_Center(axis, face);
		if((*bc).Is_Psi_N(axis, face)) return false;
		T phi=(*levelset).Phi(pos);
		return (phi>-(T)narrow_band_cell_num*mac_grid->grid.dx&&phi<(T)narrow_band_cell_num* mac_grid->grid.dx);
	}

	////Physical interface functions that defines the problem
	virtual T Off_Diag_Term(const T dt, const VectorDi& fluid_cell, const int& nbidx) const;
	virtual T Diag_Face_Term(const T dt, const int& axis, const VectorDi& face) const;
	virtual T Velocity_Correction(const T dt, const int& axis, const VectorDi& face) const;
	//Is_Valid_Cell: same as Projection<d>
	virtual bool Is_Fluid_Cell(const VectorDi& cell) const {return mac_grid->grid.Valid_Cell(cell)&&((*levelset).phi(cell)<(T)0);}
};


#endif
#endif
