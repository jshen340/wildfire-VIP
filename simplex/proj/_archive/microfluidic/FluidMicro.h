/*
Values use cgs system of units:
	 mass: g
	 velocity: cm/s
	 force: dyn (1e-5 N)
	 acceleration: Gal (1e-2 m/s^2)
	 time: s
	 length: cm
	 pressure: Ba (1e-1 Pa)
	 dynamic viscosity: P (1e-1 Pa.s)
*/
#ifndef __FluidMicro_h__
#define __FluidMicro_h__
#include "MacGrid.h"
#include "Advection.h"
#include "BoundaryCondition.h"
#include "Projection.h"
#include "KrylovSolver.h"
#include "Hashtable.h"
#include "Timer.h"
#include "Interpolation.h"
#include "FluidFunc.h"
#include "Viscosity.h"
#include "Particles.h"
#include "RandomNumber.h"

template<int d> class FluidMicro
{
	Typedef_VectorDii(d);
public:
	MacGrid<d> mac_grid;
	FaceField<real, d> velocity;
	FaceField<real, d> alpha;
	Field<ushort, d> type;		////supporting type: Fluid, Solid
	BoundaryConditionMacGrid<d> bc = BoundaryConditionMacGrid<d>(mac_grid);
	BoundaryConditionMacGridViscosity<d> bc_vis = BoundaryConditionMacGridViscosity<d>(mac_grid);
	bool verbose = false;
	/*
	Parameters of Fluid: (Gravity ignored)
	Boundary conditions set in the driver (inlet velocity, solid trap, etc)
	*/
	real rho = (real)1.0; // doesn't really matter, since our solver p*=p*dt/rho/dx
	real nu = 7e-3;	// dynamic viscosity
	real p_dia = 0.0024;	// particle diameter in cm

	TrackerPoints<d> points;
	TrackerPoints<d> p_trap;
	RandomNumber random_number;
	Array<VectorD> traps;
	Array<real> t_angs;
	Array<int> traps_count;
	Array<int> traps_open;	// open condition: 0 for open, 1 for closing, 2 for closed
	real trap_eff;

	//statistics
	real cells_sum = (real)0; 
	real sq_sum = (real)0; 
	int cell_max = INT_MIN; 
	int cell_min = INT_MAX;

	Projection<d> projection;

	FluidMicro() :projection(mac_grid, velocity, alpha, type, bc) {}

	virtual void Initialize(const VectorDi& cell_counts, const real dx, const VectorD& domain_min = VectorD::Zero())
	{
		mac_grid.Initialize(cell_counts, dx, domain_min);
		Initialize_Fields();
	}

	virtual void Initialize_Fields()
	{
		velocity.Resize(mac_grid.grid.cell_counts, (real)0);
		alpha.Resize(mac_grid.grid.cell_counts, (real)1);
		type.Resize(mac_grid.grid.cell_counts, (ushort)CellType::Fluid);
	}

	virtual void Advance(const real dt)
	{
		Timer<real> timer; timer.Reset();
		Apply_Viscosity_And_Advection(dt);
		if(verbose) timer.Elapse_And_Output_And_Reset("Advection");
		Enforce_Incompressibility();
		if(verbose) timer.Elapse_And_Output_And_Reset("Projection");
		//std::cout << "****************************************" << std::endl;
		// display number of trapped particles
		//std::cout << "Number of trapped particles: " << p_trap.Size() << std::endl;
		// display number of free particles
		//std::cout << "Number of free particles: " << points.Size() - p_trap.Size() << std::endl;
		//std::cout << "****************************************" << std::endl;
		//std::cout << "Detailed Statistics: " << std::endl;
		cells_sum = (real)0; 
		sq_sum = (real)0; 
		cell_max = INT_MIN; 
		cell_min = INT_MAX;
		trap_eff = (real)0;
		for (unsigned int i = 0; i < traps.size(); i++) {
			//std::cout << "Trap ID: " << i << " captures " << traps_count[i] << " cells." << std::endl;
			if (traps_count[i] > 0) {
				trap_eff += (real)1;
			}
			cells_sum += (real)traps_count[i];
			sq_sum += (real)traps_count[i] * (real)traps_count[i];
			if (cell_max < traps_count[i]) cell_max = traps_count[i];
			if (cell_min > traps_count[i]) cell_min = traps_count[i];
		}
		trap_eff /= traps.size();
		//std::cout << "Max Cells Captures in a Trap: " << cell_max << std::endl;
		//std::cout << "Min Cells Captures in a Trap: " << cell_min << std::endl;
		//std::cout << "Average Cells Captures in a Trap: " << cells_sum/traps.size() << std::endl;
		//std::cout << "Standard Deviation: " << sqrt(sq_sum/traps.size()- cells_sum / traps.size() * cells_sum / traps.size()) << std::endl;
		//std::cout << "****************************************" << std::endl;
	}

	virtual void Apply_Viscosity_And_Advection(const real dt)
	{
		// transfer velocity from grid to particles
		Interpolation<d> intp(mac_grid.grid);
		int p_size = points.Size();
#pragma omp parallel for
		for (int i = 0; i < p_size; i++) {
			if (!Valid_Particle(i))continue;
			VectorD v = intp.Interpolate_Face_Vectors_Cubic(velocity, points.X(i));
			points.V(i) = v;
		}

#pragma omp parallel for
		for (int i = 0; i < p_size; i++) {
			if (!Valid_Particle(i))continue;
			points.X(i) += points.V(i) * dt;
		}
		////Two cases for removing particles:
		////1. outside grid; 2. inside solid;
		////If 2, add particle to p_trap
		for (int i = 0; i < p_size; i++) {
			if (!Valid_Particle(i))continue;
			VectorD pos = points.X(i);
			VectorDi cell = mac_grid.grid.Cell_Coord(pos);
			bool outside = !mac_grid.grid.Valid_Cell(cell);

			bool trapped = false;
			for (unsigned int i = 0; i < traps.size(); i++) {
				if (Trapped(pos, traps[i])) {
					trapped = true;
					traps_count[i] += 1;
					//if (traps_count[i] >= 30 && traps_open[i] == 0) {	// change the number of saturation particles here
					//	traps_open[i] = 1;
					//}
					break;
				}
			}
			if (outside || trapped) {
				Remove_Particle(i);
				if (trapped) {
					int idx = p_trap.Add_Element();
					p_trap.X(idx) = pos;
					p_trap.V(idx) = VectorD::Zero();
				}
			}
		}
		
		FaceField<real, d> vel_copy = velocity;
		Viscosity<d> viscosity(mac_grid, velocity, alpha, type, bc_vis);
		viscosity.Solve(dt, nu);
		Advection::Semi_Lagrangian(dt, mac_grid, velocity, mac_grid, vel_copy);
		velocity = vel_copy;
		Update_Cell_Types();
	}

	virtual void Update_Cell_Types()
	{
		// update traps open/close depends on particles captured
		for (unsigned int i = 0; i < traps.size(); i++) {
			if (traps_open[i] == 1) {
				Close_Trap(traps[i]);
				traps_open[i] == 2;
			}
		}

		int cell_num = mac_grid.grid.cell_counts.prod();
#pragma omp parallel for
		for (int i = 0; i < cell_num; i++) {
			VectorDi cell = mac_grid.grid.Cell_Coord(i);
			type(cell) = bc.Is_Psi_D(cell) ? (ushort)CellType::Solid : (ushort)CellType::Fluid;
		}
	}

	virtual void Enforce_Incompressibility()
	{
		Enforce_Boundary_Conditions();
		projection.Project();
	}

	virtual void Enforce_Boundary_Conditions()
	{
		for (auto p : bc.psi_N_values) {
			int axis = p.first[0]; int face_index = p.first[1]; real value = p.second;
			velocity.face_fields[axis].array[face_index] = value;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	////helper functions

	////parallel seeding particles
	Heap<int, std::greater<int> > invalid_point_heap;
	inline bool Valid_Particle(const int i) const { return points.I(i) != -1; }
	inline int Append_Particle() { points.Add_Element(); return points.Size() - 1; }
	int Add_Particle() { if (!invalid_point_heap.empty()) { int p = invalid_point_heap.top(); invalid_point_heap.pop(); points.I(p) = 0; return p; } else return Append_Particle(); }
	void Remove_Particle(const int i) { invalid_point_heap.push(i); points.I(i) = -1; points.V(i) = VectorD::Zero(); points.X(i) = VectorD::Zero(); }

	inline void Calc_Seed_Particle_Index(const int n, Array<int>& idx)
	{
		idx.resize(n); for (int i = 0; i < n; i++) { idx[i] = Add_Particle(); }
	}

	void Seed_Particles(const VectorDi& cell, const int n)
	{
		Array<int> idx;
		Calc_Seed_Particle_Index(n, idx);
		Seed_Particles(cell, idx);
	}

	void Seed_Particles(const VectorDi& cell, const Array<int>& idx)
	{
		VectorD node_pos = mac_grid.grid.Node(cell);
		int n = (int)idx.size();
		for (int i = 0; i < n; i++) {
			int p = idx[i];
			VectorD pos = node_pos + random_number.VectorValue<d>() * mac_grid.grid.dx;
			points.X(p) = pos;
			points.V(p) = VectorD::Zero();
		}	////velocity will be interpolated from grid
	}

	// check if a particles is trapped
	bool Trapped(const VectorD& pos, const VectorD& c)
	{
		// check if in the upper half
		VectorD c_upper = c; c_upper[1] += 0.001;
		if ((pos - c_upper).norm() >= (0.0107 - p_dia / 2.0) && (pos - c_upper).norm() <= 0.0107 + p_dia / 2.0)
		{
			if (pos[0] >= c_upper[0] && pos[1] >= c[1]) return true;
		}
		// check if in the lower half
		VectorD c_lower = c; c_lower[1] -= 0.001;
		if ((pos - c_lower).norm() >= (0.0107 - p_dia / 2.0) && (pos - c_lower).norm() <= 0.0107 + p_dia / 2.0)
		{
			if (pos[0] >= c_lower[0] && pos[1] <= c[1]) return true;
		}

		return false;
	}

	void Close_Trap(const VectorD& trap)
	{
		iterate_cell(iter, mac_grid.grid) {
			VectorDi cell = iter.Coord(); VectorD pos = mac_grid.grid.Center(cell);
			bool in_trap = false;
			// check if in the upper half
			VectorD c_upper = trap; c_upper[1] += 0.001;
			if ((pos - c_upper).norm() <= 0.014)
			{
				if (pos[0] >= c_upper[0] && pos[1] >= c_upper[1]) in_trap = true;
			}
			// check if in the lower half
			VectorD c_lower = trap; c_lower[1] -= 0.001;
			if ((pos - c_lower).norm() <= 0.014)
			{
				if (pos[0] >= c_lower[0] && pos[1] <= c_lower[1]) in_trap = true;
			}
			// check if in the center gap
			if ((pos[1] <= c_upper[1] && pos[1] >= c_lower[1]) && (pos[0] >= c_lower[0] && (pos[0] - c_lower[0]) <= 0.014)) in_trap = true;
			if (in_trap) {
				bc.Set_Psi_D(cell, (ushort)CellType::Solid);
				for (int axis = 0; axis < d; axis++)for (int side = 0; side < 2; side++) { VectorDi face = mac_grid.Cell_Incident_Face(axis, cell, side); bc.Set_Psi_N(axis, face, (real)0); }
			}
		}
	}

	inline bool Is_Fluid_Cell(const VectorDi& cell) const { return type(cell) == (ushort)CellType::Fluid; }
	inline bool Is_Solid_Cell(const VectorDi& cell) const { return type(cell) == (ushort)CellType::Solid; }

	real Max_Abs(const FaceField<real, d>& field_q) const
	{
		real max_abs = (real)0;
		for (int i = 0; i < d; i++) {
			int n = (int)field_q.face_fields[i].array.size();
			for (int j = 0; j < n; j++) {
				real abs_v = abs(field_q.face_fields[i].array[j]);
				if (abs_v > max_abs)max_abs = abs_v;
			}
		}
		return max_abs;
	}
};
#endif
