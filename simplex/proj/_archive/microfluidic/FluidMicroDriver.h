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
		fluid.Advance(dt);
	}

	virtual void Write_Output_Files(const int frame)
	{
		Base::Write_Output_Files(frame);
		if (frame == 0) {
			std::string file_name = frame_dir + "/grid";
			fluid.mac_grid.grid.Write_To_File_3d(file_name);
			std::cout << "Write to file " << file_name << std::endl;
		}

		if (frame >= 1100) {
			////Write velocity
			{std::string file_name = frame_dir + "/velocity";
			fluid.velocity.Write_To_File_3d(file_name); }

			////Write BC
			{std::string file_name = frame_dir + "/psi_D";
			Particles<d> particles;
			for (auto p : fluid.bc.psi_D_values) {
				VectorDi cell = fluid.mac_grid.grid.Cell_Coord(p.first);
				VectorD pos = fluid.mac_grid.grid.Center(cell);
				int i = particles.Add_Element(); particles.X(i) = pos;
			}
			particles.Write_To_File_3d(file_name); }
			{std::string file_name = frame_dir + "/psi_N";
			Particles<d> particles;
			for (auto p : fluid.bc.psi_N_values) {
				int axis = p.first[0];
				VectorDi face = fluid.mac_grid.Face_Coord(axis, p.first[1]);
				VectorD pos = fluid.mac_grid.Face_Center(axis, face);
				int i = particles.Add_Element(); particles.X(i) = pos;
			}
			File::Write_Binary_To_File(file_name, particles);
			particles.Write_To_File_3d(file_name); }

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
			{std::string file_name = frame_dir + "/tracker_points";
			fluid.points.Write_To_File_3d_Fast(file_name); }
			{std::string file_name = frame_dir + "/trapped_points";
			fluid.p_trap.Write_To_File_3d_Fast(file_name); }
		}

		if (frame == 1200) {
			std::ofstream statistics;
			statistics.open("statistics.txt");
			statistics << "Number of trapped particles: " << fluid.p_trap.Size() << "\n";
			statistics << "Percentage of cell loss: " << (real)(fluid.points.Size() - fluid.p_trap.Size())/(real)fluid.points.Size() << "\n";
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
		case 1: {	////test
			fluid.Initialize(cell_counts, (real)length / cell_counts[0]);
			fluid.velocity.Fill((real).2, 0);

			////Source
			real source_speed = (real)2; int axis = 0;
			iterate_face_in_one_dim(axis, iter, fluid.mac_grid) {
				const VectorDi& face = iter.Coord();
				if (face[0] == 0/*||face[0]==fluid.mac_grid.face_grids[0].node_counts[0]-1*/)fluid.bc.Set_Psi_N(0, face, source_speed);
			}
			////Wall
			{axis = 1; iterate_face_in_one_dim(axis, iter, fluid.mac_grid) {
				const VectorDi& face = iter.Coord();
				if (face[axis] == 0 || face[axis] == fluid.mac_grid.face_grids[axis].node_counts[axis] - 1)fluid.bc.Set_Psi_N(axis, face, 0);
			}}
			if (d == 3) {
				axis = 2; iterate_face_in_one_dim(axis, iter, fluid.mac_grid) {
					const VectorDi& face = iter.Coord();
					if (face[axis] == 0 || face[axis] == fluid.mac_grid.face_grids[axis].node_counts[axis] - 1)fluid.bc.Set_Psi_N(axis, face, 0);
				}
			}
			////Solid
			VectorD box_center = (real).5 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); box_center[0] = (real).2;
			//VectorD box_size=VectorD::Ones()*(real).025;Box<d> obstacle(box_center-box_size,box_center+box_size);//
			real r = (real).04; Sphere<d> obstacle(box_center, r);
			iterate_cell(iter, fluid.mac_grid.grid) {
				const VectorDi& cell = iter.Coord(); const VectorD& pos = fluid.mac_grid.grid.Center(cell);
				if (obstacle.Inside(pos)) {
					fluid.bc.Set_Psi_D(cell, (ushort)CellType::Solid);
					for (int axis = 0; axis < d; axis++)for (int side = 0; side < 2; side++) { VectorDi face = fluid.mac_grid.Cell_Incident_Face(axis, cell, side); fluid.bc.Set_Psi_N(axis, face, (real)0); }
				}
			}

			fluid.projection.solver_mode = SolverType::MULTIGRID_AUTO;
			fluid.projection.update_A = false;
		}break;
		case 2: {	////microfluidic chip
			// inlet velocity: 1ul/min~=16.67ucm^3/s
			// domain: 0.24 x 0.12 (cm), res 256 x 128, left bottom corner: (0,0)
			// thickness: 0.006 (cm)
			// inlet and outlet width: 0.1/0.03 (cm) => inlet velocity=16.67e-6cm^3/s/0.006/0.1~=0.03cm/s
			// trap: inner radius: 0.0107; outer radius: 0.014; gap: 0.002 (cm)
			srand(std::time(NULL));

			length = 0.40; cell_counts[0] = s; cell_counts[1] = s / 2;
			fluid.Initialize(cell_counts, (real)length / cell_counts[0]);
			fluid.velocity.Fill((real).0, 0);

			////Source (height = 0.01 ~ 0.11 at left inlet)
			real inlet_cells = 0.1 / fluid.mac_grid.grid.dx;	// inlet height: 0.11-0.01 = 0.1, for 256x128, ~0.1/0.12*128~=107
			real source_speed = (real).03; int axis = 0;
			iterate_face_in_one_dim(axis, iter, fluid.mac_grid) {
				const VectorDi& face = iter.Coord();
				VectorD pos = fluid.mac_grid.Face_Center(axis, face);
				if (face[0] == 0 && pos[1] > 0.05 && pos[1] < 0.15) {
					fluid.bc.Set_Psi_N(axis, face, source_speed);
					VectorDi cell = MacGrid<d>::Face_Right_Cell(axis, face);
					fluid.Seed_Particles(cell, (int)pow(4, d));	// adjust to change the number of particles, take 20 traps for instance, =20000/800*100/inlet_cells~=5
				}
			}
			////Wall
			// upper and lower walls of the domain
			{axis = 1; iterate_face_in_one_dim(axis, iter, fluid.mac_grid) {
				const VectorDi& face = iter.Coord();
				if (face[axis] == 0 || face[axis] == fluid.mac_grid.face_grids[axis].node_counts[axis] - 1)fluid.bc.Set_Psi_N(axis, face, 0);
			}}
			// inner wall
			iterate_cell(iter, fluid.mac_grid.grid) {
				const VectorDi& cell = iter.Coord(); const VectorD& pos = fluid.mac_grid.grid.Center(cell);
				// bottom left solid
				if (pos[1] <= 0.05) {
					real left_boundary = 0.04 - 0.4 * pos[1];
					if (pos[0] <= left_boundary) {
						fluid.bc.Set_Psi_D(cell, (ushort)CellType::Solid);
						for (int axis = 0; axis < d; axis++)for (int side = 0; side < 2; side++) { VectorDi face = fluid.mac_grid.Cell_Incident_Face(axis, cell, side); fluid.bc.Set_Psi_N(axis, face, (real)0); }
					}
				}
				// bottom right solid
				if (pos[1] <= 0.075) {
					real right_boundary = 0.25 + 1.0 / 0.75 * pos[1];
					if (pos[0] >= right_boundary) {
						fluid.bc.Set_Psi_D(cell, (ushort)CellType::Solid);
						for (int axis = 0; axis < d; axis++)for (int side = 0; side < 2; side++) { VectorDi face = fluid.mac_grid.Cell_Incident_Face(axis, cell, side); fluid.bc.Set_Psi_N(axis, face, (real)0); }
					}
				}
				// upper left solid
				if (pos[1] >= 0.15) {
					real left_boundary = 0.02 + 0.4 * (pos[1] - 0.15);
					if (pos[0] <= left_boundary) {
						fluid.bc.Set_Psi_D(cell, (ushort)CellType::Solid);
						for (int axis = 0; axis < d; axis++)for (int side = 0; side < 2; side++) { VectorDi face = fluid.mac_grid.Cell_Incident_Face(axis, cell, side); fluid.bc.Set_Psi_N(axis, face, (real)0); }
					}
				}
				// upper right solid
				if (pos[1] >= 0.125) {
					real right_boundary = 0.35 - 1.0/0.75 * (pos[1]-0.125);
					if (pos[0] >= right_boundary) {
						fluid.bc.Set_Psi_D(cell, (ushort)CellType::Solid);
						for (int axis = 0; axis < d; axis++)for (int side = 0; side < 2; side++) { VectorDi face = fluid.mac_grid.Cell_Incident_Face(axis, cell, side); fluid.bc.Set_Psi_N(axis, face, (real)0); }
					}
				}
			}
			////Traps
			// manually design
			//column 1
			VectorD trap_c1 = (real).90 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c1[0] = 0.044; trap_c1[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c1[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_1 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c2 = (real).65 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c2[0] = 0.044; trap_c2[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c2[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_2 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c3 = (real).40 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c3[0] = 0.044; trap_c3[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c3[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_3 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c4 = (real).10 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c4[0] = 0.044; trap_c4[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c4[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_4 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;

			//column 2
			VectorD trap_c5 = (real).77 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c5[0] = 0.094; trap_c5[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c5[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_5 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c6 = (real).52 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c6[0] = 0.094; trap_c6[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c6[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_6 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c7 = (real).25 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c7[0] = 0.094; trap_c7[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c7[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_7 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			
			// column 3
			VectorD trap_c8 = (real).90 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c8[0] = 0.134; trap_c8[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c8[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_8 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c9 = (real).65 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c9[0] = 0.134; trap_c9[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c9[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_9 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c10 = (real).40 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c10[0] = 0.134; trap_c10[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c10[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_10 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c11 = (real).10 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c11[0] = 0.134; trap_c11[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c11[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_11 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			
			//column 4
			VectorD trap_c12 = (real).90 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c12[0] = 0.164; trap_c12[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c12[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_12 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c13 = (real).65 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c13[0] = 0.164; trap_c13[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c13[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_13 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c14 = (real).40 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c14[0] = 0.164; trap_c14[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c14[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_14 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c15 = (real).10 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c15[0] = 0.164; trap_c15[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c15[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_15 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;

			//column 5
			VectorD trap_c16 = (real).72 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c16[0] = 0.194; trap_c16[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c16[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_16 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c17 = (real).47 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c17[0] = 0.194; trap_c17[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c17[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_17 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c18 = (real).25 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c18[0] = 0.194; trap_c18[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c18[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_18 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;

			//Column 6
			VectorD trap_c19 = (real).90 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c19[0] = 0.214; trap_c19[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c19[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_19 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c20 = (real).65 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c20[0] = 0.214; trap_c20[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c20[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_20 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c21 = (real).40 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c21[0] = 0.214; trap_c21[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c21[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_21 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c22 = (real).10 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c22[0] = 0.214; trap_c22[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c22[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_22 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;

			//column 7
			VectorD trap_c23 = (real).75 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c23[0] = 0.234; trap_c23[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c23[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_23 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c24 = (real).55 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c24[0] = 0.234; trap_c24[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c24[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_24 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c25 = (real).35 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c25[0] = 0.234; trap_c25[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c25[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_25 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c26 = (real).15 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c26[0] = 0.234; trap_c26[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c26[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_26 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;

			//column 8
			//VectorD trap_c27 = (real).87 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c27[0] = 0.264;
			VectorD trap_c28 = (real).65 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c28[0] = 0.254; trap_c28[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c28[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_28 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c29 = (real).45 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c29[0] = 0.254; trap_c29[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c29[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_29 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c30 = (real).20 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c30[0] = 0.254; trap_c30[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c30[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_30 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;

			//column 9
			VectorD trap_c31 = (real).71 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c31[0] = 0.274; trap_c31[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c31[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_31 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c32 = (real).49 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c32[0] = 0.274; trap_c32[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c32[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_32 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c33 = (real).27 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c33[0] = 0.274; trap_c33[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c33[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_33 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;

			// column 10
			VectorD trap_c34 = (real).68 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c34[0] = 0.294; trap_c34[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c34[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_34 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c35 = (real).32 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c35[0] = 0.294; trap_c35[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c35[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_35 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;

			// column 11
			VectorD trap_c37 = (real).63 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c37[0] = 0.314; trap_c37[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c37[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_37 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			VectorD trap_c38 = (real).42 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c38[0] = 0.314; trap_c38[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c38[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_38 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;

			//column 12
			VectorD trap_c39 = (real).49 * (fluid.mac_grid.grid.domain_min + fluid.mac_grid.grid.domain_max); trap_c39[0] = 0.334; trap_c39[0] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; trap_c39[1] += ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * fluid.mac_grid.grid.dx * (real)2; real ang_39 = ((rand() / double(RAND_MAX)) * 2.0 - 1.0) * atan(1) * 4.0 / 9.0;
			
			
			fluid.traps.push_back(trap_c1);  fluid.t_angs.push_back(ang_1); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c2);  fluid.t_angs.push_back(ang_2); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c3);  fluid.t_angs.push_back(ang_3); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c4);  fluid.t_angs.push_back(ang_4); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c5);  fluid.t_angs.push_back(ang_5); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c6);  fluid.t_angs.push_back(ang_6); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c7);  fluid.t_angs.push_back(ang_7); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c8);  fluid.t_angs.push_back(ang_8); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c9);  fluid.t_angs.push_back(ang_9); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c10);  fluid.t_angs.push_back(ang_10); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c11);  fluid.t_angs.push_back(ang_11); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c12);  fluid.t_angs.push_back(ang_12); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c13);  fluid.t_angs.push_back(ang_13); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c14);  fluid.t_angs.push_back(ang_14); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c15);  fluid.t_angs.push_back(ang_15); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c16);  fluid.t_angs.push_back(ang_16); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c17);  fluid.t_angs.push_back(ang_17); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c18);  fluid.t_angs.push_back(ang_18); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c19);  fluid.t_angs.push_back(ang_19); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c20);  fluid.t_angs.push_back(ang_20); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c21);  fluid.t_angs.push_back(ang_21); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c22);  fluid.t_angs.push_back(ang_22); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c23);  fluid.t_angs.push_back(ang_23); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c24);  fluid.t_angs.push_back(ang_24); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c25);  fluid.t_angs.push_back(ang_25); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c26);  fluid.t_angs.push_back(ang_26); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			//fluid.traps.push_back(trap_c27);  fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c28);  fluid.t_angs.push_back(ang_28); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c29);  fluid.t_angs.push_back(ang_29); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c30);  fluid.t_angs.push_back(ang_30); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c31);  fluid.t_angs.push_back(ang_31); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c32);  fluid.t_angs.push_back(ang_32); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c33);  fluid.t_angs.push_back(ang_33); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c34);  fluid.t_angs.push_back(ang_34); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c35);  fluid.t_angs.push_back(ang_35); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			//fluid.traps.push_back(trap_c36);  fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c37);  fluid.t_angs.push_back(ang_37); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c38);  fluid.t_angs.push_back(ang_38); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			fluid.traps.push_back(trap_c39);  fluid.t_angs.push_back(ang_39); fluid.traps_count.push_back(0); fluid.traps_open.push_back(0);
			iterate_cell(iter, fluid.mac_grid.grid) {
				VectorDi cell = iter.Coord(); VectorD pos = fluid.mac_grid.grid.Center(cell);
				for (unsigned int i = 0; i < fluid.traps.size(); i++) {
					if (In_Trap(pos, fluid.traps[i], fluid.t_angs[i])) {
						fluid.bc.Set_Psi_D(cell, (ushort)CellType::Solid);
						for (int axis = 0; axis < d; axis++)for (int side = 0; side < 2; side++) { VectorDi face = fluid.mac_grid.Cell_Incident_Face(axis, cell, side); fluid.bc.Set_Psi_N(axis, face, (real)0); }
					}
				}
			}

			fluid.projection.solver_mode = SolverType::MULTIGRID_AUTO;
			fluid.projection.update_A = false;

			// write traps to disk
			std::ofstream trap_info;
			trap_info.open("trap_info.txt");
			for (unsigned int i = 0; i < fluid.traps.size(); i++) {
				trap_info << "trap: " << i << ": [" << fluid.traps[i][0] << "," << fluid.traps[i][1] << "], " << fluid.t_angs[i] << "\n";
			}
			trap_info.close();
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
};

#endif
