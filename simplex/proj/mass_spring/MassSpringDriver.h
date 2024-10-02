//#####################################################################
// Mass spring driver
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com, Yankai Mao
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __MassSpringDriver_h__
#define __MassSpringDriver_h__
#include <memory>
#include "Common.h"
#include "Driver.h"
#include "Mesh.h"
#include "MeshFunc.h"
#include "SoftBodyMassSpring.h"

template<int d> class MassSpringDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;using Base::verbose;
public:
	SoftBodyMassSpring<d> soft_body;
	std::shared_ptr<SegmentMesh<d> > segment_mesh=nullptr;
	std::shared_ptr<TriangleMesh<d> > triangle_mesh=nullptr;
	std::shared_ptr<TetrahedronMesh<d> > tetrahedron_mesh = nullptr;

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		soft_body.Advance(dt);
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);
		
		////Write Particles
		{std::string file_name=frame_dir+"/particles";
		soft_body.particles.Write_To_File_3d(file_name);}

		////Write springs
		if(segment_mesh!=nullptr){
			std::string file_name=frame_dir+"/segment_mesh";
			segment_mesh->Write_To_File_3d(file_name);}

		if(triangle_mesh!=nullptr){
			std::string file_name=frame_dir+"/triangle_mesh";
			triangle_mesh->Write_To_File_3d(file_name);}

		if (tetrahedron_mesh != nullptr) {
			std::string file_name = frame_dir + "/tetrahedron_mesh";
			tetrahedron_mesh->Write_To_File_3d(file_name);
		}

		////Write BC
		{std::string file_name=frame_dir+"/psi_D";
		Particles<d> psi_D_particles;
		for(auto p:soft_body.psi_D_values){int idx=p.first;
			int i=psi_D_particles.Add_Element();psi_D_particles.X(i)=soft_body.particles.X(idx);}
		psi_D_particles.Write_To_File_3d(file_name);}
	}
	
	virtual void Initialize()
	{
		verbose=true;

		switch(test){
		case 1:{	////rope with one fixed point
			real length=(real)1;int n=20;real step=length/(real)n;
			soft_body.particles.Resize(n);
			segment_mesh.reset(new SegmentMesh<d>(soft_body.particles.XPtr()));

			for(int i=0;i<n;i++){
				soft_body.particles.X(i)=VectorD::Unit(0)*(real)i*step;
				soft_body.particles.M(i)=(real)1;}
			for(int i=0;i<n-1;i++){Vector2i s(i,i+1);
				segment_mesh->Elements().push_back(s);}

			soft_body.Set_Psi_D(0);

			Array<Vector2i> edges;MeshFunc::Get_Edges(*segment_mesh,edges);
			soft_body.springs=edges;
		}break;
		case 2:{	////cloth with two fixed points
			real length=(real)1;int w=3*scale;int h=3*scale;real step=length/(real)w;
			TriangleMesh<d> tri_mesh_copy;
			MeshFunc::Initialize_Herring_Bone_Mesh(w,h,step,&tri_mesh_copy,0,2);
			int n=(int)tri_mesh_copy.Vertices().size();
			
			soft_body.particles.Resize(n);
			triangle_mesh.reset(new TriangleMesh<d>(soft_body.particles.XPtr()));	
			triangle_mesh->elements=tri_mesh_copy.elements;

			for(int i=0;i<n;i++){
				soft_body.particles.X(i)=tri_mesh_copy.Vertices()[i];
				soft_body.particles.M(i)=(real)1;}

			soft_body.Set_Psi_D(0);
			soft_body.Set_Psi_D(w-1);

			Array<Vector2i> edges;MeshFunc::Get_Edges(*triangle_mesh,edges);
			soft_body.springs=edges;
		}break;
		case 3: {	////box (can not be compiled when d equals to 2)
			TetrahedronMesh<d> tet_mesh_copy;
			VectorDi cell_counts=scale*VectorDi::Ones();real dx=(real)1/(real)scale;
			cell_counts[1]/=4;cell_counts[2]/=4;
			MeshFunc::Initialize_Lattice_Mesh(cell_counts,dx,&tet_mesh_copy);
			int n=(int)tet_mesh_copy.Vertices().size();

			soft_body.particles.Resize(n);
			tetrahedron_mesh.reset(new TetrahedronMesh<d>(soft_body.particles.XPtr()));
			tetrahedron_mesh->elements=tet_mesh_copy.elements;

			for(int i=0;i<n;i++){
				soft_body.particles.X(i)=tet_mesh_copy.Vertices()[i];
				soft_body.particles.M(i)=(real)1;}

			int axis=0;
			for(int i=0;i<(int)tet_mesh_copy.Vertices().size();i++){
				const VectorD& pos=tet_mesh_copy.Vertices()[i];
				if(pos[axis]==0){soft_body.Set_Psi_D(i);;}}

			Array<Vector2i> edges;MeshFunc::Get_Edges(*tetrahedron_mesh, edges);
			soft_body.springs=edges;
		}break;
		case 4:{	////cloth with 4 fixed points
			real length=(real)1;int w=20*scale;int h=20*scale;real step=length/(real)w;
			TriangleMesh<d> tri_mesh_copy;
			MeshFunc::Initialize_Herring_Bone_Mesh(w,h,step,&tri_mesh_copy,0,2);
			int n=(int)tri_mesh_copy.Vertices().size();
			
			soft_body.particles.Resize(n);
			triangle_mesh.reset(new TriangleMesh<d>(soft_body.particles.XPtr()));	
			triangle_mesh->elements=tri_mesh_copy.elements;

			for(int i=0;i<n;i++){
				soft_body.particles.X(i)=tri_mesh_copy.Vertices()[i];
				soft_body.particles.M(i)=(real)1;}

			soft_body.Set_Psi_D(0);
			soft_body.Set_Psi_D(w-1);
			soft_body.Set_Psi_D((h-1) * w);
			soft_body.Set_Psi_D(h * w - 1);

			Array<Vector2i> edges;MeshFunc::Get_Edges(*triangle_mesh,edges);
			soft_body.springs=edges;
		}break;
		case 5:{	////rope with 2 fixed points
			real length=(real)1;int n=20;real step=length/(real)n;
			soft_body.particles.Resize(n);
			segment_mesh.reset(new SegmentMesh<d>(soft_body.particles.XPtr()));

			for(int i=0;i<n;i++){
				soft_body.particles.X(i)=VectorD::Unit(0)*(real)i*step;
				soft_body.particles.M(i)=(real)1;}
			for(int i=0;i<n-1;i++){Vector2i s(i,i+1);
				segment_mesh->Elements().push_back(s);}

			soft_body.Set_Psi_D(0);
			soft_body.Set_Psi_D(n-1);

			Array<Vector2i> edges;MeshFunc::Get_Edges(*segment_mesh,edges);
			soft_body.springs=edges;
		}break;
		case 6:{	////sloping rod using bending spring
			real length=(real)1;int n=30;real step=length/(real)n;
			soft_body.particles.Resize(n);
			segment_mesh.reset(new SegmentMesh<d>(soft_body.particles.XPtr()));

			VectorD dir((real)0.5,(real)0.866,(real).0);

			for(int i=0;i<n;i++){
				soft_body.particles.X(i)=dir*(real)i*step;
				soft_body.particles.M(i)=(real)1 / (i+1);}
			
			for (int i = 0; i < n*7/10;i++)
				soft_body.Set_Psi_D(i);

			for(int i=0;i<n-1;i++){Vector2i s(i,i+1);
				segment_mesh->Elements().push_back(s);}

			Array<Vector2i> edges;MeshFunc::Get_Edges(*segment_mesh,edges);
			soft_body.springs=edges;
			soft_body.use_bending_spring_line=true;

			soft_body.use_newmark=false;
		}break;
		case 7:{	////stif cloth with fixed edge using bending spring
			real length=(real)1;int w=10*scale;int h=10*scale;real step=length/(real)w;
			TriangleMesh<d> tri_mesh_copy;
			MeshFunc::Initialize_Herring_Bone_Mesh(w,h,step,&tri_mesh_copy,0,2);
			int n=(int)tri_mesh_copy.Vertices().size();
			
			soft_body.particles.Resize(n);
			triangle_mesh.reset(new TriangleMesh<d>(soft_body.particles.XPtr()));	
			triangle_mesh->elements=tri_mesh_copy.elements;

			for(int i=0;i<n;i++){
				soft_body.particles.X(i)=tri_mesh_copy.Vertices()[i];
				soft_body.particles.M(i)=(real)1/w/h;}

			for (int i = 0; i < w; i++)
				for (int j = 0;j < 2;j++) {
					soft_body.Set_Psi_D(w*j+i);
				}

			Array<Vector2i> edges;MeshFunc::Get_Edges(*triangle_mesh,edges);
			soft_body.springs=edges;
			soft_body.use_bending_spring_surface=true;
			soft_body.use_newmark=false;
		}break;
		case 9:{	//// very stiff cloth using dihedral constraint
			real length=(real)1;int w=20*scale;int h=20*scale;real step=length/(real)w;
			TriangleMesh<d> tri_mesh_copy;
			MeshFunc::Initialize_Herring_Bone_Mesh(w,h,step,&tri_mesh_copy,0,2);
			int n=(int)tri_mesh_copy.Vertices().size();
			
			soft_body.particles.Resize(n);
			triangle_mesh.reset(new TriangleMesh<d>(soft_body.particles.XPtr()));	
			triangle_mesh->elements=tri_mesh_copy.elements;

			for(int i=0;i<n;i++){
				soft_body.particles.X(i)=tri_mesh_copy.Vertices()[i];
				soft_body.particles.M(i)=(real)1/w/h;}

			for (int i = 0; i < w; i++)
				for (int j = 0;j < 2;j++) {
					soft_body.Set_Psi_D(w*j+i);
				}

			Array<Vector2i> edges;MeshFunc::Get_Edges(*triangle_mesh,edges);
			soft_body.springs=edges;
			soft_body.use_bending_force=true;
			soft_body.use_newmark=true;
			soft_body.ks_0=1e3;
			soft_body.kd_0=1e1;
		}break;
		}

		soft_body.Initialize();
	}
};

#endif