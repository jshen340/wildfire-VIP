//#####################################################################
// Rigid body driver
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __RigidDriver_h__
#define __RigidDriver_h__
#include <fstream>
#include "Common.h"
#include "File.h"
#include "Particles.h"
#include "Driver.h"
#include "Rigid.h"
#include "SpatialHashing.h"
#include "Mesh.h"
#include "LevelSet.h"
#include "RandomNumber.h"

template<int d> class RigidDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;
public:
	Array<Rigid<d> > rigids;
	Particles<d> particles;
	Array<Vector2i> particle_indices;
	VectorD g=VectorD::Unit(1)*(real)-1;
	real dx=(real).15;
	
	SpatialHashing<d,16> spatial_hashing;
	Array<Vector2i> particle_particle_collision_pairs;
	Array<Vector2i> particle_environment_collision_pairs;
	bool use_env_collision=false;
	std::set<int> fixed_rigids;
	std::map<int,int> vis_mesh_vtx;
	
	Array<SegmentMesh<d> > obj_mesh;
	SegmentMesh<d> vis_mesh;

	Array<LevelSet<d> > levelsets;
	RandomNumber random=RandomNumber((real)-1,(real)1);

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		for(int i=0;i<(int)rigids.size();i++){
			rigids[i].Clear_Force_And_Torque();
		
			////add random forces
			for(int j=particle_indices[i][0];j<particle_indices[i][1]-1;j++){
				VectorD f=random.VectorValue<d>();
				rigids[i].Add_Force_In_World_Space(f,particles.X(j));}}

		//Collision_Detection();
		//Collision_Response();

		for(int i=0;i<(int)rigids.size();i++){
			if(fixed_rigids.find(i)!=fixed_rigids.end())continue;
			rigids[i].Advance(dt);
			rigids[i].Sync_Particles(particles,particle_indices[i][0],particle_indices[i][1]);
		}
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);

		//Write particles
		{std::string file_name=frame_dir+"/particles";
		particles.Write_To_File_3d(file_name);}

		//Write mesh
		Sync_Object_Mesh_To_Vis_Mesh();
		{std::string file_name=frame_dir+"/segment_mesh";
		vis_mesh.Write_To_File_3d(file_name);}

		std::cout<<"write "<<frame<<std::endl;
	}
	
	void Add_Body_Force(const int idx,const VectorD& g)
	{
		for(int i=particle_indices[idx][0];i<particle_indices[idx][1];i++){
			rigids[idx].Add_Force_In_World_Space(particles.M(i)*g,particles.X(i));}
	}

	void Collision_Detection()
	{
		spatial_hashing.Clear_Voxels();
		spatial_hashing.Update_Voxels(particles.XRef());

		particle_environment_collision_pairs.clear();
		if(use_env_collision)for(int idx=0;idx<(int)rigids.size();idx++)
			for(int i=particle_indices[idx][0];i<particle_indices[idx][1];i++){
				if(particles.X(i)[1]<(real)0){
					particle_environment_collision_pairs.push_back(Vector2i(i,0));}}

		particle_particle_collision_pairs.clear();
		for(int i=0;i<particles.Size();i++){
			ArrayF<int,64> nbs;spatial_hashing.Find_Nbs(particles.X(i),particles.XRef(),dx,nbs);
			for(int k=1;k<=nbs[0];k++){int j=nbs[k];
				if(particles.I(i)==particles.I(j)||i>=j)continue;
				if((particles.X(i)-particles.X(j)).norm()<dx){
					particle_particle_collision_pairs.push_back(Vector2i(i,j));}}}
	}

	void Collision_Response()
	{
		real ks=(real)1e4;
		real kd=(real)1e2;

		////collision with ground
		if(use_env_collision)for(auto& cp:particle_environment_collision_pairs){
			int i=cp[0];int r_i=particles.I(i);
			VectorD dir=VectorD::Unit(1);
			real length=particles.X(i)[1];
			VectorD rel_dir=VectorD::Zero()-particles.V(i);
			VectorD f_s=ks*(length)*(-dir);
			VectorD f_d=kd*rel_dir.dot(dir)*dir;
			rigids[r_i].Add_Force_In_World_Space(f_s+f_d,particles.X(i));}

		////collision among objects
		for(auto& cp:particle_particle_collision_pairs){
			int i=cp[0];int j=cp[1];int r_i=particles.I(i);int r_j=particles.I(j);
			real rest_length=dx;
			VectorD dir=particles.X(j)-particles.X(i);
			real length=dir.norm();
			if(length>rest_length)continue;
			dir/=length;
			VectorD rel_dir=particles.V(j)-particles.V(i);
			VectorD f_s=ks*(length-rest_length)*dir;
			VectorD f_d=kd*rel_dir.dot(dir)*dir;
			rigids[r_i].Add_Force_In_World_Space(f_s+f_d,particles.X(i));
			rigids[r_j].Add_Force_In_World_Space(-f_s-f_d,particles.X(j));}
	}

	virtual void Initialize()
	{
		switch(test){
		case 1:{	////three rigid body falling
			int w=8;int h=4;dx=(real)1/(real)w;int obj_p_num=w*h;int obj_n=3;int p_n=obj_p_num*obj_n;
			particles.Resize(p_n);
			rigids.resize(obj_n);
			particle_indices.resize(obj_n);
			
			for(int i=0;i<obj_n;i++){
				VectorD pos=VectorD::Unit(1)*(real)(i+1)+VectorD::Unit(0)*.5*(real)(i);
				for(int ii=0;ii<w;ii++){
					for(int jj=0;jj<h;jj++){
						int idx=i*obj_p_num+ii*h+jj;
						VectorD p=pos+VectorD::Unit(0)*dx*(real)ii+VectorD::Unit(1)*dx*(real)jj;
						particles.X(idx)=p;}}
				for(int j=0;j<obj_p_num;j++){
					particles.M(i*obj_p_num+j)=(real)1;
					particles.F(i*obj_p_num+j)=VectorD::Zero();
					particles.I(i*obj_p_num+j)=i;
					particles.V(i*obj_p_num+j)=VectorD::Zero();}}

			for(int i=0;i<obj_n;i++){
				Vector2i indices=Vector2i(i*obj_p_num,i*obj_p_num+obj_p_num);
				rigids[i].Initialize_From_Particles(particles,indices[0],indices[1]);
				particle_indices[i]=indices;}
		}break;
		case 2:{	////a rigid body from mesh
			Initialize_Objects_From_Files();
		}break;
		case 3:{	////a 3D cube rigid body 
			if constexpr (d==3){
				VectorDi counts=Vector3i::Ones()*4;real dx=(real)1/(real)counts[0];
				Grid<d> grid(counts,dx);
				int pn=grid.cell_counts.prod();
				real density=(real)1;
				real m0=density*pow(dx,d);
				particles.Resize(pn);
				iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();
					int i=grid.Cell_Index(cell);
					particles.M(i)=m0;
					particles.F(i)=VectorD::Zero();
					particles.I(i)=i;
					particles.V(i)=VectorD::Zero();
					particles.X(i)=grid.Center(cell);}
			
				rigids.resize(1);
				rigids[0].Initialize_From_Particles(particles,0,pn);
				
				particle_indices.resize(1);
				particle_indices[0]=Vector2i(0,pn);
			}
		}break;
		}

		spatial_hashing.Initialize(dx);
		cfl=(real).001;
	}

	void Initialize_Objects_From_Files()
	{
		//std::string file_name="two-link-chain-vertices.txt";
		//std::string file_name="50-link-open-chain.txt";
		std::string file_name="30x30-grid.txt";

		std::ifstream in(file_name);

		////read each obj mesh from file
		int vtx_n;
		while(in>>vtx_n){
			std::cout<<"vtx_n: "<<vtx_n<<std::endl;
			SegmentMesh<d> mesh;
			for(int i=0;i<vtx_n;i++){
				VectorD v;in>>v[0]>>v[1];
				(*mesh.vertices).push_back(v);}
			for(int i=0;i<vtx_n-1;i++){
				mesh.elements.push_back(Vector2i(i,i+1));}
			mesh.elements.push_back(Vector2i(vtx_n-1,0));
			obj_mesh.push_back(mesh);}

		////create simulation particles
		real dxx=dx*(real).5;
		particle_indices.resize(obj_mesh.size());
		rigids.resize(obj_mesh.size());
		int mesh_vtx_n=0;
		vis_mesh_vtx.clear();
		for(int r_i=0;r_i<obj_mesh.size();r_i++){
			int pidx_start=particles.Size();
			SegmentMesh<d>& mesh=obj_mesh[r_i];
			int rigid_body_p_n=0;
			for(int i=0;i<mesh.elements.size();i++){
				Vector2i e=mesh.elements[i];
				VectorD p0=(*mesh.vertices)[e[0]];
				VectorD p1=(*mesh.vertices)[e[1]];
				real length=(p1-p0).norm();
				VectorD dir=(p1-p0)/length;
				int n=(int)(length/dxx);
				real step=length/(real)n;
				int start=particles.Add_Elements(n);
				vis_mesh_vtx[e[0]+mesh_vtx_n]=rigid_body_p_n;
				rigid_body_p_n+=n;
				for(int j=0;j<n;j++){
					particles.X(start+j)=p0+dir*(real)j*step;
					particles.M(start+j)=(real)1;
					particles.F(start+j)=VectorD::Zero();
					particles.I(start+j)=r_i;}}
			int pidx_end=particles.Size();
			particle_indices[r_i]=Vector2i(pidx_start,pidx_end);
			rigids[r_i].Initialize_From_Particles(particles,pidx_start,pidx_end);
			mesh_vtx_n+=(int)(*mesh.vertices).size();
		}

		////sync to vis mesh
		for(int r_i=0;r_i<obj_mesh.size();r_i++){
			int n=(int)(*vis_mesh.vertices).size();
			SegmentMesh<d>& mesh=obj_mesh[r_i];
			for(int j=0;j<(*mesh.vertices).size();j++){
				(*vis_mesh.vertices).push_back((*mesh.vertices)[j]);}
			for(int j=0;j<mesh.elements.size();j++){
				Vector2i e=mesh.elements[j];e[0]+=n;e[1]+=n;
				vis_mesh.elements.push_back(e);}}
	}

	void Sync_Object_Mesh_To_Vis_Mesh()
	{
		int idx=0;
		for(int r_i=0;r_i<obj_mesh.size();r_i++){
			for(int j=0;j<(*obj_mesh[r_i].vertices).size();j++){int p_idx=vis_mesh_vtx[idx];
				(*vis_mesh.vertices)[idx]=rigids[r_i].World_Space_Point(rigids[r_i].r0[p_idx]);idx++;}}
	}

	void Initialize_Levelsets()
	{
		levelsets.resize(rigids.size());
		for(int i=0;i<rigids.size();i++){
			Box<d> box=MeshFunc::Bounding_Box(rigids.r0);
			VectorDi cell_counts=VectorDi::Ones()*32;
			VectorD domain_min=box.min_corner;dx=box.Edge_Lengths().Max()/(real)32;
			Grid<d> grid=Grid<d>(cell_counts,dx,domain_min);
			levelsets[i].Initialize(grid);}
	}
};

#endif