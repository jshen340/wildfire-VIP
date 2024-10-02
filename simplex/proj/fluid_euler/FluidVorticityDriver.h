//////////////////////////////////////////////////////////////////////////
// Euler fluid vorticity driver
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __FluidVorticityDriver_h__
#define __FluidVorticityDriver_h__
#include "Common.h"
#include "File.h"
#include "GeometryPrimitives.h"
#include "FluidEulerOld.h"
#include "Particles.h"
#include "MarchingCubes.h"
#include "Driver.h"
#include "Vorticity.h"
#include "FluidFunc.h"
#include "AnalyticalFields.h"

template<int d> class FluidVorticityDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;
public:
	FluidEulerOld<d> fluid;

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		//fluid.Advance(dt);
		//Check_Vorticity_Conversion();
	}

	void Test_Vorticity_To_Velocity()
	{
		MacGrid<d>& mac_grid=fluid.mac_grid;
		FaceField<real,d>& velocity=fluid.velocity;velocity.Resize(fluid.mac_grid.grid.cell_counts);
		Field<real,d> vorticity_on_nodes(fluid.mac_grid.grid.node_counts,(real)0);

		iterate_node(iter,fluid.mac_grid.grid){const VectorDi& node=iter.Coord();
			const VectorD pos=fluid.mac_grid.grid.Node(node);
			vorticity_on_nodes(node)=AnaField::Taylor_Green_Vorticity(pos);}

		Vorticity::Vorticity_To_Velocity(mac_grid,vorticity_on_nodes,velocity);
	}

	void Check_Vorticity_Conversion()
	{
		MacGrid<d>& mac_grid=fluid.mac_grid;
		FaceField<real,d>& velocity=fluid.velocity;

		FaceField<real,d> u0(fluid.mac_grid.grid.cell_counts);
		iterate_face(axis,iter,fluid.mac_grid){const VectorDi& face=iter.Coord();
			const VectorD pos=fluid.mac_grid.Face_Center(axis,face);
			VectorD v=AnaField::Taylor_Green_Velocity(pos);
			u0(axis,face)=v[axis];}

		Field<real,d> vorticity(fluid.mac_grid.grid.node_counts,(real)0);
		FluidFunc::Curl_On_Node(fluid.mac_grid,u0,vorticity);

		iterate_node(iter,fluid.mac_grid.grid){const VectorDi& node=iter.Coord();
			const VectorD pos=fluid.mac_grid.grid.Node(node);
			real w=AnaField::Taylor_Green_Vorticity(pos);
			std::cout<<"("<<w<<", "<<vorticity(node)<<")";
			vorticity(node)=w;
			//if(node[0]==1||node[1]==1||
			//	node[0]==fluid.mac_grid.grid.node_counts[0]-2||
			//	node[1]==fluid.mac_grid.grid.node_counts[1]-2)
			//	std::cout<<"vor "<<node.transpose()<<": "<<vorticity(node)<<std::endl;
		}

		velocity.Resize(fluid.mac_grid.grid.cell_counts);
		//velocity=u0;

		////Vorticity::Vorticity_To_Velocity(fluid.mac_grid,vorticity,velocity);

		Poisson<d> poisson;
		Grid<d> node_grid(mac_grid.grid.node_counts,mac_grid.grid.dx,mac_grid.grid.domain_min-VectorD::Ones()*mac_grid.grid.dx*(real).5);
		poisson.Initialize(node_grid);
		poisson.use_multigrid_solver=false;
		iterate_cell(iter,node_grid){const VectorDi& cell=iter.Coord();
			poisson.rhs(cell)=-vorticity(cell);}

		////boundary conditions
		iterate_face(axis,iter,poisson.mac_grid){const VectorDi& face=iter.Coord();
			if(/*face[axis]==0||*/face[axis]==poisson.mac_grid.face_grids[axis].node_counts[axis]-1){
				VectorD face_pos=poisson.mac_grid.Face_Center(axis,face);
				VectorD v=AnaField::Taylor_Green_Velocity(face_pos);
				if(axis==0)poisson.Set_Psi_N(axis,face,-v[1]);
				else poisson.Set_Psi_N(axis,face,v[0]);}}

		poisson.Build_And_Solve();
		FluidFunc::Curl_On_Face(mac_grid,poisson.p,velocity);
		
		FluidFunc::Curl_On_Node(fluid.mac_grid,velocity,vorticity);
		iterate_node(iter,fluid.mac_grid.grid){const VectorDi& node=iter.Coord();
			const VectorD pos=fluid.mac_grid.grid.Node(node);
			if(node[0]==1||node[1]==1||
				node[0]==fluid.mac_grid.grid.node_counts[0]-2||
				node[1]==fluid.mac_grid.grid.node_counts[1]-2)
				std::cout<<"vorxx "<<node.transpose()<<": "<<vorticity(node)<<std::endl;
		}

		std::cout<<"poisson.nc: "<<poisson.nc<<std::endl;

		iterate_face(axis,iter,fluid.mac_grid){const VectorDi& face=iter.Coord();
			std::cout<<"["<<u0(axis,face)<<", "<<velocity(axis,face)<<"] ";}
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);
		if(frame==0){
			std::string file_name=frame_dir+"/grid";
			fluid.mac_grid.grid.Write_To_File_3d(file_name);
			std::cout<<"Write to file "<<file_name<<std::endl;}
		
		////Write velocity
		{std::string file_name=frame_dir+"/velocity";
		fluid.velocity.Write_To_File_3d(file_name);}

		////Write BC
		{std::string file_name=frame_dir+"/psi_D";
		Particles<d> particles;
		for(auto p:fluid.bc.psi_D_values){
			VectorDi cell=fluid.mac_grid.grid.Cell_Coord(p.first);
			VectorD pos=fluid.mac_grid.grid.Center(cell);
			int i=particles.Add_Element();particles.X(i)=pos;}
		particles.Write_To_File_3d(file_name);}
		{std::string file_name=frame_dir+"/psi_N";
		Particles<d> particles;
		for(auto p:fluid.bc.psi_N_values){int axis=p.first[0];
			VectorDi face=fluid.mac_grid.Face_Coord(axis,p.first[1]);
			VectorD pos=fluid.mac_grid.Face_Center(axis,face);
			int i=particles.Add_Element();particles.X(i)=pos;}
		File::Write_Binary_To_File(file_name,particles);
		particles.Write_To_File_3d(file_name);}

		////Write fluid type
		{std::string file_name=frame_dir+"/fluid_type";
		Field<real,d> fluid_type;fluid_type.Resize(fluid.mac_grid.grid.cell_counts);
		iterate_cell(iter,fluid.mac_grid.grid){const VectorDi& cell=iter.Coord();
			if(fluid.type(cell)==(ushort)CellType::Fluid)fluid_type(cell)=(real)0;
			else fluid_type(cell)=(real)1;}
		fluid_type.Write_To_File_3d(file_name);}
	}
	
	virtual void Initialize()
	{
		int s=scale;real length=pi;VectorDi cell_counts=VectorDi::Ones()*s;
		
		fluid.Initialize(cell_counts,(real)length/cell_counts[0]);
		//Check_Vorticity_Conversion();
		Test_Vorticity_To_Velocity();

		//cell_counts[1]=cell_counts[0]/2;if(d>2)cell_counts[2]=cell_counts[0]/2;

		//switch(test){
		//case 1:{
		//	fluid.Initialize(cell_counts,(real)length/cell_counts[0]);
		//	fluid.velocity.Fill((real).2,0);

		//	////Source
		//	real source_speed=(real)1;int axis=0;
		//	iterate_face_in_one_dim(axis,iter,fluid.mac_grid){const VectorDi& face=iter.Coord();
		//		if(face[0]==0||face[0]==fluid.mac_grid.face_grids[0].node_counts[0]-1)fluid.bc.Set_Psi_N(0,face,source_speed);}
		//	//////Wall
		//	//axis=1;iterate_face_in_one_dim(axis,iter,fluid.mac_grid){const VectorDi& face=iter.Coord();
		//	//	if(face[axis]==0||face[axis]==fluid.mac_grid.face_grids[axis].node_counts[axis]-1)fluid.bc.Set_Psi_N(axis,face,0);}
		//	//////Solid
		//	//VectorD box_center=(real).5*(fluid.mac_grid.grid.domain_min+fluid.mac_grid.grid.domain_max);
		//	//VectorD box_size=VectorD::Ones()*(real).1;Box<d> box(box_center-box_size,box_center+box_size);
		//	//iterate_cell(iter,fluid.mac_grid.grid){const VectorDi& cell=iter.Coord();const VectorD& pos=fluid.mac_grid.grid.Center(cell);
		//	//	if(box.Inside(pos)){fluid.bc.Set_Psi_D(cell,(ushort)CellType::Solid);
		//	//		for(int axis=0;axis<d;axis++)for(int side=0;side<2;side++){VectorDi face=fluid.mac_grid.Cell_Incident_Face(axis,cell,side);fluid.bc.Set_Psi_N(axis,face,(real)0);}}}
		//}break;
		//case 2:{
		//	fluid.Initialize(cell_counts,(real)length/cell_counts[0]);
		//	fluid.velocity.Fill((real).2,0);

		//	real source_speed=(real)1;
		//	VectorD box_center=fluid.mac_grid.grid.Center();
		//	box_center[0]=(real).2;VectorD box_size=VectorD::Ones()*(real).1;box_size[0]=fluid.mac_grid.grid.dx*(real)2;
		//	Box<d> box(box_center-box_size,box_center+box_size);

		//	////Initialize source
		//	VectorD source_velocity=VectorD::Unit(0)*(real)1;
		//	iterate_cell(iter,fluid.mac_grid.grid){const VectorDi& cell=iter.Coord();
		//		const VectorD& pos=fluid.mac_grid.grid.Center(cell);
		//			if(box.Inside(pos)){fluid.bc.Set_Psi_D(cell,(ushort)CellType::Solid);
		//			for(int axis=0;axis<d;axis++)for(int side=0;side<2;side++){
		//				VectorDi face=fluid.mac_grid.Cell_Incident_Face(axis,cell,side);fluid.bc.Set_Psi_N(axis,face,source_velocity[axis]);}}}
		//}break;}
		//fluid.Enforce_Boundary_Conditions();
	
	}
};

#endif