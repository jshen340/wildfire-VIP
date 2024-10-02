//#####################################################################
// Cloth driver
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __ClothDriver_h__
#define __ClothDriver_h__
#include <memory>
#include "Common.h"
#include "Driver.h"
#include "Mesh.h"
#include "Cloth.h"

class ClothDriver : public Driver
{Typedef_VectorDii(3);using Base=Driver;using Base::verbose;
public:
	ClothMassSpring cloth;
	std::shared_ptr<SegmentMesh<3> > segment_mesh=nullptr;

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		cloth.Advance(dt);
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);
		
		////Write Particles
		{std::string file_name=frame_dir+"/particles";
		cloth.particles.Write_To_File_3d(file_name);}

		////Write springs
		if(segment_mesh!=nullptr){
			std::string file_name=frame_dir+"/segment_mesh";
			segment_mesh->Write_To_File_3d(file_name);}

		////Write BC
		{std::string file_name=frame_dir+"/psi_D";
		Particles<3> psi_D_particles;
		for(auto p:cloth.psi_D_values){int idx=p.first;
			int i=psi_D_particles.Add_Element();psi_D_particles.X(i)=cloth.particles.X(idx);}
		psi_D_particles.Write_To_File_3d(file_name);}
	}
	
	virtual void Initialize()
	{
		verbose=true;

		switch(test){
		case 1:{	////Cloth
			real length=(real)1;int w=8*scale;int h=12*scale;real dx=length/(real)w;
			Grid<2> grid(Vector2i(w,h),dx);
			cloth.grid=grid;
			int n=grid.node_counts.prod();
			cloth.particles.Resize(n);

			iterate_node_d(iter,grid,2){
				Vector2i node=iter.Coord();int node_idx=grid.Node_Index(node);Vector2 pos=grid.Node(node);
				cloth.particles.X(node_idx)=Vector3(pos[0],(real)0,pos[1]);
				cloth.particles.M(node_idx)=(real)1;}

			cloth.Set_Psi_D(grid.Node_Index(Vector2i(0,0)));
			cloth.Set_Psi_D(grid.Node_Index(Vector2i(w,0)));

			cloth.Initialize_Springs();
		}break;
		}

		cloth.Initialize();

		////Initialize segement_mesh for visualization only
		segment_mesh.reset(new SegmentMesh<3>(cloth.particles.XPtr()));
		segment_mesh->elements=cloth.springs;
	}
};

#endif