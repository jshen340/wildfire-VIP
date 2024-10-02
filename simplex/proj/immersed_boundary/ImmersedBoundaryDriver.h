//////////////////////////////////////////////////////////////////////////
// Mesh To Level Set Driver
// Copyright (c) (2021-), Fan Feng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __ImmersedBoundaryDriver_h__
#define __ImmersedBoundaryDriver_h__
#include "Common.h"
#include "Driver.h"
#include "Field.h"
#include "Grid.h"
#include "File.h"
#include "Mesh.h"
#include "ImmersedBoundary.h"
#include "MeshAdvFunc.h"
#include "MeshFunc.h"
#include "MarchingCubes.h"
#include <iostream>
#include <fstream>
#ifdef USE_TINY_OBJ_LOADER
#include "TinyObjLoader.h"
#endif

template<int d> class ImmersedBoundaryDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;
public:
	ImmersedBoundary<d> fluid;
	
	virtual real CFL() const
	{
		real epsilon = (real)1e-5;
		return cfl * fluid.fluid.mac_grid.grid.dx / (fluid.fluid.Max_Abs(fluid.fluid.velocity) + epsilon);
	}

	virtual void Initialize()
	{
		VectorDi cell_counts=VectorDi::Ones()*scale;
		real dx=(real)1/(real)scale;
		fluid.Initialize(cell_counts, dx);
		
		switch(test){
			case 1:{
				real step_size = dx / (real)2;
				auto& soft_body = fluid.soft_body;

				if constexpr (d == 2) {
					real m0 = (real)1.; // Cannot be less than one
					int n = 16; real length = n * step_size;
					soft_body.particles.Resize(n);
					std::shared_ptr<SegmentMesh<d> > segment_mesh = nullptr;
					segment_mesh.reset(new SegmentMesh<d>(soft_body.particles.XPtr()));

					VectorD start_point = fluid.fluid.mac_grid.grid.Position(AuxFunc::V<d>((real).5 - length * (real).5, (real).95, (real).5));
					for (int i = 0; i < n; i++) {
						soft_body.particles.X(i) = start_point + VectorD::Unit(0) * (real)i * step_size;
						soft_body.particles.M(i) = m0;
						soft_body.particles.F(i) = VectorD::Zero();
					}
					for (int i = 0; i < n - 1; i++) {
						Vector2i s(i, i + 1);
						segment_mesh->Elements().push_back(s);
					}

					////bending
					for (int i = 0; i < n - 2; i++) {
						Vector2i s(i, i + 1);
						segment_mesh->Elements().push_back(s);
					}
					for (int i = 0; i < n - 3; i++) {
						Vector2i s(i, i + 2);
						segment_mesh->Elements().push_back(s);
					}

					Array<Vector2i> edges; MeshFunc::Get_Edges(*segment_mesh, edges);
					soft_body.springs = edges;
				}
				else if constexpr (d == 3) {
					real m0 = (real)1.;
					int w = 16; int h = 16; real length = w * step_size;

					std::shared_ptr<TriangleMesh<d> > triangle_mesh = nullptr;

					TriangleMesh<d> tri_mesh_copy;
					MeshFunc::Initialize_Herring_Bone_Mesh(w, h, step_size, &tri_mesh_copy, 0, 2);
					int n = (int)tri_mesh_copy.Vertices().size();

					soft_body.particles.Resize(n);
					triangle_mesh.reset(new TriangleMesh<d>(soft_body.particles.XPtr()));
					triangle_mesh->elements = tri_mesh_copy.elements;

					VectorD start_point = fluid.fluid.mac_grid.grid.Position(AuxFunc::V<d>((real).5 - length * (real).5, (real).95, (real).5 - length * (real).5));

					for (int i = 0; i < n; i++) {
						soft_body.particles.X(i) = start_point + tri_mesh_copy.Vertices()[i];
						soft_body.particles.M(i) = m0;
					}

					Array<Vector2i> edges; MeshFunc::Get_Edges(*triangle_mesh, edges);
					soft_body.springs = edges;
				}
				
				soft_body.ks_0 = (real)1e4;
				soft_body.kd_0 = (real)1e2;
				fluid.soft_body.g = VectorD::Unit(1) * (real)-9.8;

				fluid.soft_body.Initialize();

				//set walls
				for (int axis = 0; axis < d; axis++)
					iterate_face_in_one_dim(axis, iter, fluid.fluid.mac_grid) {
						const VectorDi& face = iter.Coord();
						if (face[axis] == 0 || face[axis] == fluid.fluid.mac_grid.face_grids[axis].node_counts[axis] - 1){
							fluid.fluid.bc.Set_Psi_N(axis, face, 0); 
						}
				}
			}break;
		}
	}

	virtual void Advance_One_Time_Step(const real dt, const real time) {
		fluid.Advance(dt);
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);
		if (frame == 0) {
			std::string file_name = frame_dir + "/grid";
			fluid.fluid.mac_grid.grid.Write_To_File_3d(file_name);
			std::cout << "Write to file " << file_name << std::endl;
		}

		////Write velocity
		{std::string file_name = frame_dir + "/velocity";
		fluid.fluid.velocity.Write_To_File_3d(file_name); }

		{std::string file_name = frame_dir + "/solid_particles";
		fluid.soft_body.particles.Write_To_File_3d(file_name);}

		if (frame == 0)
		{
			std::string file_name = frame_dir + "/psi_N";
			Particles<d> particles;
			for (auto p : fluid.fluid.bc.psi_N_values) {
				int axis = p.first[0];
				VectorDi face = fluid.fluid.mac_grid.Face_Coord(axis, p.first[1]);
				VectorD pos = fluid.fluid.mac_grid.Face_Center(axis, face);
				int i = particles.Add_Element(); particles.X(i) = pos;
			}
			File::Write_Binary_To_File(file_name, particles);
			particles.Write_To_File_3d(file_name);
		}

		std::cout << "Write to frame " << frame << std::endl;
	}
};
#endif