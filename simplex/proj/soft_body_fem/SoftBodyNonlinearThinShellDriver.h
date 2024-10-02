//#####################################################################
// Soft Body Nonlinear Thin Shell Driver
// Copyright (c) (2021-), Fan Feng, fan.feng.gr@dartmouth.edu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __SoftBodyNonlinearThinShellDynamicDriver_h__
#define __SoftBodyNonlinearThinShellDynamicDriver_h__
#include <fstream>
#include "Common.h"
#include "Field.h"
#include "Driver.h"
#include "AuxFunc.h"
#include "MeshFunc.h"
#include "SoftBodyNonlinearFemThinShell.h"

template<int d> class SoftBodyNonlinearThinShellDriver : public Driver
{Typedef_VectorDii(d);Typedef_MatrixD(d);using Base=Driver;
public:
	SurfaceMesh<d> mesh;
	SoftBodyNonlinearFemThinShell<d> thin_shell;
	bool use_quasi_static=false;

	virtual void Initialize()
	{
		Initialize_Mesh();
		Initialize_Boundary_Conditions();
		Initialize_Materials();
	}

	void Initialize_Mesh()
	{
		if constexpr (d == 2) {
			switch (test) {
			case 1: {
				MeshFunc::Initialize_Segment_Mesh<d>(VectorD::Zero(), VectorD::Unit(0), scale, &mesh,false,true);
			}break;
			}
		}

		if constexpr (d == 3) {
			switch(test){
			case 1: {
				real length = (real)1; int w = scale; real step = length / (real)w;
				MeshFunc::Initialize_Herring_Bone_Mesh(w, w, step, &mesh, 0, 2);
			}break;
			case 2: {
				real length = (real)1; int w = scale; real step = length / (real)w;
				MeshFunc::Initialize_Herring_Bone_Mesh(w, w, step, &mesh, 0, 2);
			}break;
			case 3: {
				real length = (real)1; int w = scale; real step = length / (real)w;
				MeshFunc::Initialize_Herring_Bone_Mesh(w, w, step, &mesh, 0, 2);
			}break;
			case 4: {
				real length = (real)1; int w = scale; real step = length / (real)w;
				MeshFunc::Initialize_Herring_Bone_Mesh(w, w, step, &mesh, 0, 2);
			}break;
			case 5: {
				real length = (real)1; int w = scale; real step = length / (real)w;
				MeshFunc::Initialize_Herring_Bone_Mesh(w, w, step, &mesh, 0, 2);
			}break;
			case 6:
			case 7:
			case 8: {
				real length = (real)1; int w = scale; real step = length / (real)w;
				MeshFunc::Initialize_Herring_Bone_Mesh(w, w * 2 - 1, step, &mesh, 0, 2);
				for (int i = w; i < (w * 2 - 2) * w; i++) {
					mesh.Vertices()[i] += -(real)step * 0.1 * (w - 1 - abs((int)(i / w) - (w - 1))) * VectorD::Unit(1);
				}
				//for (int i = w * (w - 1); i < w * w; i++){mesh.Vertices()[i] += -(real)0.1*VectorD::Unit(1); }
			}break;
			case 9: {
				real length = (real)1; int w = scale; real step = length / (real)w;
				MeshFunc::Initialize_Herring_Bone_Mesh(w, w * 2 - 1, step, &mesh, 0, 2);

				/*for (int i = w; i < (w * 2 - 2) * w; i++) {
					mesh.Vertices()[i] += -(real)step * 0.1 * (w - 1 - abs((int)(i / w) - (w - 1))) * VectorD::Unit(1);
				}*/
				for (int i = w * (w - 1); i < w * w; i++) { mesh.Vertices()[i] += (real)step * VectorD::Unit(1); }
				for (int i = w * (w - 3); i < w * (w - 2); i++) { mesh.Vertices()[i] += -(real)step * VectorD::Unit(1); }
				for (int i = w * (w + 1); i < w * (w + 2); i++) { mesh.Vertices()[i] += -(real)step * VectorD::Unit(1); }
			}break;
			}
		}

		thin_shell.Initialize(mesh);
	}

	void Initialize_Boundary_Conditions()
	{
		switch(test){
		case 1:{	////beam with one end fixed under gravity
			if constexpr (d == 2) {
				thin_shell.Set_Fixed(0);
			}
			if constexpr (d == 3) {
				real length = (real)1; int w = scale; real step = length / (real)w;

				for (int i = 0; i < w; i++) {
					thin_shell.Set_Fixed(i);
				}
			}
			thin_shell.use_body_force=true;
		}break;
		case 2:{	////beam with two ends fixed under gravity
			real length = (real)1; int w = scale; real step = length / (real)w;
			for (int i = 0; i < w; i++) {thin_shell.Set_Fixed(i);}
			for (int i = thin_shell.particles.Size()-w; i < thin_shell.particles.Size(); i++) {thin_shell.Set_Fixed(i);}

			thin_shell.use_body_force=true;
		}break;
		case 3:{	////pull one side and push on the other side
			real length = (real)1; int w = scale; real step = length / (real)w; real strength = 0.1;
			for (int i = 0; i < w; i++) {thin_shell.Set_Force(i, strength*VectorD::Unit(1)/(real)w);}
			for (int i = thin_shell.particles.Size()-w; i < thin_shell.particles.Size(); i++) {thin_shell.Set_Force(i,-strength * VectorD::Unit(1)/(real)w);}

			thin_shell.use_body_force=false;
		}break;
		case 4:{	////Pull two ends out
			real length = (real)1; int w = scale; real step = length / (real)w; real strength = 10;
			for (int i = 0; i < w; i++) {thin_shell.Set_Force(i,-strength*VectorD::Unit(2) / (real)w);}
			for (int i = thin_shell.particles.Size()-w; i < thin_shell.particles.Size(); i++) {thin_shell.Set_Force(i, strength * VectorD::Unit(2) / (real)w);}

			thin_shell.use_body_force=false;
		}break;
		case 5: {  ////one diagonal push and the other diagonal pull, not supported in static solve
			int w = scale;
			real strength = 0.01;
			thin_shell.Set_Force(0, strength*VectorD::Unit(1));
			thin_shell.Set_Force(thin_shell.Vtx_Num()-1, strength * VectorD::Unit(1));
			thin_shell.Set_Force(w-1, -strength * VectorD::Unit(1));
			thin_shell.Set_Force(thin_shell.Vtx_Num()-w, -strength * VectorD::Unit(1));
			thin_shell.use_body_force = false;
		}break;
		case 6: {  ////fix one side and pull the other side on the plane
			int w = scale;
			real strength = 10;
			for (int i = 0; i < w; i++) {thin_shell.Set_Fixed(i);}
			for (int i = thin_shell.particles.Size() - w; i < thin_shell.particles.Size(); i++) { thin_shell.Set_Force(i,strength*VectorD::Unit(2)/(real)w); }
			thin_shell.use_body_force = false;
		}break;
		case 7: {  ////fix one side and one corner
			int w = scale;
			real strength = 0.1;
			for (int i = 0; i < 2*w; i++) { thin_shell.Set_Fixed(i); }
			for (int i = thin_shell.particles.Size() - w; i < thin_shell.particles.Size(); i++) { thin_shell.Set_Force(i, strength * -VectorD::Unit(1) / (real)w); }
			thin_shell.use_body_force = false;
		}break;
		case 8: {  ////push one sides inward on the plane, with a crease in the middle
			int w = scale;
			real strength = 0.012;
			for (int i = 0; i < 2*w; i++) { thin_shell.Set_Fixed(i); }
			for (int i = thin_shell.particles.Size() - w; i < thin_shell.particles.Size(); i++) { thin_shell.Set_Force(i, strength * -VectorD::Unit(2)); }
			thin_shell.use_body_force = false;
		}break;
		case 9: {  ////push two sides inward on the plane, with three crease in the middle
			int w = scale;
			real strength = 0.012;
			for (int i = 0; i < 2*w; i++) { thin_shell.Set_Fixed(i); }
			for (int i = thin_shell.particles.Size() - w; i < thin_shell.particles.Size(); i++) { thin_shell.Set_Displacement(i, -0.1*VectorD::Unit(2)); }
			thin_shell.use_body_force = false;
		}break;
		}
	}

	void Initialize_Materials()
	{
		switch (test) {
			case 1:
			case 2:
			case 3:
			case 4:
			case 5:
			case 6:
			case 7:{
				thin_shell.materials.clear();
				thin_shell.Add_Material((real)1e2, (real).35);
				AuxFunc::Fill(thin_shell.material_id, 0);
				thin_shell.Initialize_Material();
			}break;
			case 8: {
				thin_shell.materials.clear();
				thin_shell.Add_Material((real)100, (real).35);
				AuxFunc::Fill(thin_shell.material_id, 0);
				std::fill(thin_shell.hs.begin(), thin_shell.hs.end(), (real)1);
				int w = scale;
				for (int i = w * (w - 1); i < w * w; i++) { thin_shell.hs[i] = (real)0.05; }
				thin_shell.Initialize_Material();
			}break;
			case 9: {
				thin_shell.materials.clear();
				thin_shell.Add_Material((real)100, (real).35);
				AuxFunc::Fill(thin_shell.material_id, 0);
				std::fill(thin_shell.hs.begin(), thin_shell.hs.end(), (real)1);
				int w = scale;
				for (int i = w * (w - 1); i < w * w; i++) { thin_shell.hs[i] = (real)0.07; }
				for (int i = w * (w - 3); i < w * (w - 2); i++) { thin_shell.hs[i] = (real)0.08; }
				for (int i = w * (w + 1); i < w * (w + 2); i++) { thin_shell.hs[i] = (real)0.06; }
				thin_shell.Initialize_Material();
			}break;
		}
	}

	virtual void Advance_One_Time_Step(const real dt, const real time)
	{
		if (use_quasi_static) {
			thin_shell.Advance_Quasi_Static();
			Write_Output_Files(1);
			exit(0);
		}
		else {
			thin_shell.Advance(dt, time);
		}
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);

		{std::string file_name=frame_dir+(d==2?"/segment_mesh":"/triangle_mesh");
		thin_shell.mesh->Write_To_File_3d(file_name);}

		{std::string file_name=frame_dir+"/particles";
		thin_shell.particles.Write_To_File_3d(file_name);}

		{std::string file_name=frame_dir+"/mat";
		int n=(int)thin_shell.material_id.size();
		Field<real,1> mat;mat.Resize(n);
		for(int i=0;i<n;i++){
			mat.array[i]=(real)thin_shell.material_id[i];
		}
		mat.Write_To_File_3d(file_name);}

		{std::string file_name = frame_dir + "/psi_D";
		Particles<d> psi_D_particles;
		for (auto& p : thin_shell.bc.psi_D_values) {
			int idx = p.first;
			int i = psi_D_particles.Add_Element(); psi_D_particles.X(i) = thin_shell.particles.X(idx);
		}
		psi_D_particles.Write_To_File_3d(file_name); }

		{std::string file_name = frame_dir + "/psi_N";
		Particles<d> psi_N_particles;
		for (auto& p : thin_shell.bc.forces) {
			int idx = p.first;
			int i = psi_N_particles.Add_Element(); psi_N_particles.X(i) = thin_shell.particles.X(idx);
			psi_N_particles.F(i) = p.second;
		}
		psi_N_particles.Write_To_File_3d(file_name); }

		if (frame == last_frame) {
			std::string file_name = output_dir + "/energy.txt";
			std::ofstream energy_file(file_name);
			for (int count = 0; count < thin_shell.energies_n.size(); count++) {
				energy_file << thin_shell.energies_n[count] << "\n";
			}
		}
	}
};

#endif