#ifndef __FluidMicroDriver_h__
#define __FluidMicroDriver_h__
#include "Common.h"
#include "File.h"
#include "GeometryPrimitives.h"
#include "FluidMicro.h"
#include "Particles.h"
#include "MarchingCubes.h"
#include "Driver.h"
#include "LevelSet.h"
#include "RenderFunc.h"
#include <time.h>
#include <fstream>

template<int d> class FluidMicroDriver : public Driver
{
	Typedef_VectorDii(d); using Base = Driver;
public:
	FluidMicro<d> fluid;

	real CFL() const
	{
		real epsilon = (real)1e-5;
		return cfl * fluid.mac_grid.grid.dx / (fluid.Max_Abs(fluid.velocity) + epsilon);
	}

	virtual void Advance_One_Time_Step(const real dt, const real time)
	{
		fluid.Advance(dt, time);
	}

	virtual void Write_Output_Files(const int frame)
	{
		Base::Write_Output_Files(frame);
		if (frame == 0) {
			std::string file_name = frame_dir + "/grid";
			fluid.mac_grid.grid.Write_To_File_3d(file_name);
			std::cout << "Write to file " << file_name << std::endl;
			// write traps to disk
			std::ofstream trap_info;
			trap_info.open("trap_info.txt");
			for (unsigned int i = 0; i < fluid.traps.size(); i++) {
				trap_info << "trap: " << i << ": [" << fluid.traps[i][0] << "," << fluid.traps[i][1] << "]" << "\n";
			}
			trap_info.close();
		}

		if (frame >= 0) {
			////Write velocity
			{std::string file_name = frame_dir + "/velocity";
			fluid.velocity.Write_To_File_3d(file_name); }

			////Write BC
			RenderFunc::Write_Dirichlet_Boundary_Conditions(frame_dir + "/psi_D", fluid.mac_grid, fluid.bc);
			RenderFunc::Write_Neumann_Boundary_Conditions(frame_dir + "/psi_N", fluid.mac_grid, fluid.bc);

			////Write fluid type
			{std::string file_name = frame_dir + "/fluid_type";
			Field<real, d> fluid_type; fluid_type.Resize(fluid.mac_grid.grid.cell_counts);
			iterate_cell(iter, fluid.mac_grid.grid) {
				const VectorDi& cell = iter.Coord();
				if (fluid.type(cell) == (ushort)CellType::Fluid)fluid_type(cell) = (real)0;
				else fluid_type(cell) = (real)1;
			}
			fluid_type.Write_To_File_3d(file_name); }

			////Write Particles
			RenderFunc::Write_Points_Float<d, real>(frame_dir + "/tracker_points", fluid.points.XRef());
			RenderFunc::Write_Points_Float<d, real>(frame_dir + "/trapped_points", fluid.p_trap.XRef());

			std::ofstream statistics;
			std::string file_name = frame_dir + "/statistics.txt";
			statistics.open(file_name);
			statistics << "Number of trapped particles: " << fluid.p_trap.Size() << "\n";
			statistics << "Percentage of cell loss: " << (real)(fluid.points.Size() - fluid.p_trap.Size()) / (real)fluid.points.Size() << "\n";
			statistics << "Trap efficiency: " << fluid.trap_eff << "\n";
			statistics << "Max Cells Captures in a Trap: " << fluid.cell_max << "\n";
			statistics << "Min Cells Captures in a Trap: " << fluid.cell_min << "\n";
			statistics << "Average Cells Captures in a Trap: " << fluid.cells_sum / fluid.traps.size() << "\n";
			statistics << "Standard Deviation: " << sqrt(fluid.sq_sum / fluid.traps.size() - fluid.cells_sum / fluid.traps.size() * fluid.cells_sum / fluid.traps.size()) << "\n";
			statistics.close();
		}

		std::cout << "Write to frame " << frame << std::endl;
	}

	virtual void Update_Face_Fraction(LevelSet<d>& _ls_solid)
	{
		iterate_face(axis, iter, fluid.mac_grid) {
			const VectorDi& face = iter.Coord();
			VectorDi node_0 = fluid.mac_grid.Face_Incident_Node(axis, face, 0);
			VectorDi node_1 = fluid.mac_grid.Face_Incident_Node(axis, face, 1);
			real phi_0 = _ls_solid.Phi(fluid.mac_grid.grid.Node(node_0));
			real phi_1 = _ls_solid.Phi(fluid.mac_grid.grid.Node(node_1));
			// no need to multiply dx
			if (phi_0 < (real)0 && phi_1 >(real)0) {
				fluid.alpha.face_fields[axis](face) = phi_0 / (phi_0 - phi_1);
			}
			else if (phi_0 > (real)0 && phi_1 < (real)0) {
				fluid.alpha.face_fields[axis](face) = phi_1 / (phi_1 - phi_0);
			}
			else if (phi_0 > (real)0 && phi_1 > (real)0) {
				fluid.alpha.face_fields[axis](face) = (real)0;
			}
			else {
				fluid.alpha.face_fields[axis](face) = (real)1;
			}
		}
	}

	virtual void Initialize()
	{
		int s = scale; real length = (real)1; VectorDi cell_counts = VectorDi::Ones(); cell_counts[0] = s;
		cell_counts[1] = cell_counts[0] / 2; if (d > 2)cell_counts[2] = cell_counts[0] / 2;

		switch (test) {
		case 1: {	
			srand(std::time(NULL));

			length = 0.20; cell_counts[0] = s; cell_counts[1] = s / 2;
			fluid.Initialize(cell_counts, (real)length / cell_counts[0]);
			fluid.velocity.Fill((real).0, 0);

			////Wall
			int axis = 0;
			// upper and lower walls of the domain
			{axis = 1; iterate_face_in_one_dim(axis, iter, fluid.mac_grid) {
				const VectorDi& face = iter.Coord();
				if (face[axis] == 0 || face[axis] == fluid.mac_grid.face_grids[axis].node_counts[axis] - 1)fluid.bc.Set_Psi_N(axis, face, 0);
			}}
			// set left inlet and right outlet
			{axis = 0; iterate_face_in_one_dim(axis, iter, fluid.mac_grid) {
				const VectorDi& face = iter.Coord();
				VectorD face_center = fluid.mac_grid.Face_Center(axis, face);
				if (face[axis] == 0 && (face_center[1] <= fluid.inlet_offset || face_center[1] >= (fluid.inlet_offset+fluid.inlet_width)))fluid.bc.Set_Psi_N(axis, face, 0);
				if (face[axis] == fluid.mac_grid.face_grids[axis].node_counts[axis] - 1 && (face_center[1] <= fluid.outlet_offset || face_center[1] >= (fluid.outlet_offset + fluid.outlet_width)))fluid.bc.Set_Psi_N(axis, face, 0);
			}}
			////Source 
			axis = 0;
			iterate_face_in_one_dim(axis, iter, fluid.mac_grid) {
				const VectorDi& face = iter.Coord();
				VectorD face_center = fluid.mac_grid.Face_Center(axis, face);
				if (face[axis] == 0 && (face_center[1] >= fluid.inlet_offset && face_center[1] <= (fluid.inlet_offset + fluid.inlet_width))) {
					fluid.bc.Set_Psi_N(axis, face, fluid.source_speed);
					VectorDi cell = MacGrid<d>::Face_Right_Cell(axis, face);
					fluid.Seed_Particles(cell, fluid.seed_count);
				}
			}

			VectorD c1 = VectorD::Zero(); c1[0] = 0.05; c1[1] = 0.02;
			fluid.traps.push_back(c1); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			VectorD c2 = VectorD::Zero(); c2[0] = 0.05; c2[1] = 0.05;
			fluid.traps.push_back(c2); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			VectorD c3 = VectorD::Zero(); c3[0] = 0.05; c3[1] = 0.08;
			fluid.traps.push_back(c3); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			iterate_cell(iter, fluid.mac_grid.grid) {
				VectorDi cell = iter.Coord(); VectorD pos = fluid.mac_grid.grid.Center(cell);
				for (unsigned int i = 0; i < fluid.traps.size(); i++) {
					if (In_Trap_new_1(pos, fluid.traps[i])) {
						fluid.bc.Set_Psi_D(cell, (ushort)CellType::Solid);
						for (int axis = 0; axis < d; axis++)for (int side = 0; side < 2; side++) { VectorDi face = fluid.mac_grid.Cell_Incident_Face(axis, cell, side); fluid.bc.Set_Psi_N(axis, face, (real)0); }
					}
				}
			}
			fluid.projection.update_A = true;
		}break;
		}
		fluid.Enforce_Boundary_Conditions();
	}

	// check if a grid cell is located inside a trap centered at c with angle ang
	// dimensions of the trap are defined in the comments above
	bool In_Trap(const VectorD& pos, const VectorD& c, const real& ang)
	{
		//real ang = -atan(1) * 4.0 / 6.0;
		// check if in the upper half
		VectorD c_upper = c; c_upper[0] -= 0.001 * sin(ang); c_upper[1] += 0.001 * cos(ang);
		if ((pos - c_upper).norm() >= 0.0107 && (pos - c_upper).norm() <= 0.014)
		{
			Vector2d bound1(cos(ang), sin(ang)); Vector2d bound2(-sin(ang),cos(ang));
			VectorD pos_vec = pos - c_upper;
			if ((pos_vec[0]*bound1[0]+pos_vec[1]*bound1[1]) >= 0 && (pos_vec[0] * bound2[0] + pos_vec[1] * bound2[1]) >= 0) return true;
		}
		// check if in the lower half
		VectorD c_lower = c; c_lower[0] += 0.001 * sin(ang); c_lower[1] -= 0.001 * cos(ang);
		if ((pos - c_lower).norm() >= 0.0107 && (pos - c_lower).norm() <= 0.014)
		{
			Vector2d bound1(cos(ang), sin(ang)); Vector2d bound2(sin(ang), -cos(ang));
			VectorD pos_vec = pos - c_lower;
			if ((pos_vec[0] * bound1[0] + pos_vec[1] * bound1[1]) >= 0 && (pos_vec[0] * bound2[0] + pos_vec[1] * bound2[1]) >= 0) return true;
		}
		return false;
	}

	bool In_Trap_new_1(const VectorD& pos, const VectorD& c)
	{
		// check if in the upper half
		VectorD up_min = c; up_min[1] += 0.004; VectorD up_max = c; up_max[0] += 0.012; up_max[1] += 0.006;
		if (pos[0] >= up_min[0] && pos[0] <= up_max[0] && pos[1] >= up_min[1] && pos[1] <= up_max[1]) return true;
		// check if in the lower half
		VectorD bt_min = c; bt_min[1] -= 0.006; VectorD bt_max = c; bt_max[0] += 0.012; bt_max[1] -= 0.004;
		if (pos[0] >= bt_min[0] && pos[0] <= bt_max[0] && pos[1] >= bt_min[1] && pos[1] <= bt_max[1]) return true;
		// check if in the right upper
		VectorD ru_min = c; ru_min[0] += 0.01; ru_min[1] += 0.001; VectorD ru_max = c; ru_max[0] += 0.012; ru_max[1] += 0.004;
		if (pos[0] >= ru_min[0] && pos[0] <= ru_max[0] && pos[1] >= ru_min[1] && pos[1] <= ru_max[1]) return true;
		// check if in the right lower
		VectorD rl_min = c; rl_min[0] += 0.01; rl_min[1] -= 0.004; VectorD rl_max = c; rl_max[0] += 0.012; rl_max[1] -= 0.001;
		if (pos[0] >= rl_min[0] && pos[0] <= rl_max[0] && pos[1] >= rl_min[1] && pos[1] <= rl_max[1]) return true;
		return false;
	}
	bool In_Trap_new_2(const VectorD& pos, const VectorD& c)
	{
		// check if in the upper half
		VectorD up_min = c; up_min[1] += 0.004; VectorD up_max = c; up_max[0] += 0.008; up_max[1] += 0.006;
		if (pos[0] >= up_min[0] && pos[0] <= up_max[0] && pos[1] >= up_min[1] && pos[1] <= up_max[1]) return true;
		// check if in the lower half
		VectorD bt_min = c; bt_min[1] -= 0.006; VectorD bt_max = c; bt_max[0] += 0.008; bt_max[1] -= 0.004;
		if (pos[0] >= bt_min[0] && pos[0] <= bt_max[0] && pos[1] >= bt_min[1] && pos[1] <= bt_max[1]) return true;
		// check if in the right upper
		VectorD ru_min = c; ru_min[0] += 0.01; ru_min[1] += 0.001; VectorD ru_max = c; ru_max[0] += 0.012; ru_max[1] += 0.006;
		if (pos[0] >= ru_min[0] && pos[0] <= ru_max[0] && pos[1] >= ru_min[1] && pos[1] <= ru_max[1]) return true;
		// check if in the right lower
		VectorD rl_min = c; rl_min[0] += 0.01; rl_min[1] -= 0.006; VectorD rl_max = c; rl_max[0] += 0.012; rl_max[1] -= 0.001;
		if (pos[0] >= rl_min[0] && pos[0] <= rl_max[0] && pos[1] >= rl_min[1] && pos[1] <= rl_max[1]) return true;
		return false;
	}
};

#endif
