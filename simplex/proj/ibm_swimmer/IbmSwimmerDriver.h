//#####################################################################
// Immersed boundary swimmer driver
// Author: Xingzhe He
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __IbmSwimmerDriver_h__
#define __IbmSwimmerDriver_h__
#include "Common.h"
#include "File.h"
#include "FaceField.h"
#include "GeometryPrimitives.h"
#include "FluidEuler.h"
#include "Driver.h"
#include "IbmSwimmer.h"

template<int d> class IbmSwimmerDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;
public:
	FluidEuler<d> fluid;
	SoftSwimmer<d> swimmer;
	FaceField<real,d> coupling_forces;

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		swimmer.Advance(dt,time);
		fluid.Advection(dt);
		Solid_Fluid_Coupling(dt);
		iterate_face(axis,iter,fluid.mac_grid){const VectorDi& face=iter.Coord();
			fluid.velocity(axis,face)+=coupling_forces(axis,face)*dt;}	////TODO: divide by mass
		fluid.Enforce_Incompressibility();
	}

	void Solid_Fluid_Coupling(const real dt)
	{
		coupling_forces.Fill((real)0);
		FaceField<real,d> face_weights;face_weights.Resize(fluid.mac_grid.grid.cell_counts,(real)0);
		Interpolation<d> intp(fluid.mac_grid);
		
		real coef=(real)1e-3;
		for(int i=0;i<swimmer.particles.Size();i++){
			VectorD v_s=swimmer.particles.V(i);
			VectorD pos=swimmer.particles.X(i);
			VectorD v_f=intp.Interpolate_Face_Vectors(fluid.velocity,pos);
			VectorD f=coef*swimmer.particles.M(i)*(v_s-v_f)/dt;
			intp.Interpolate_Point_To_Faces(pos,f,coupling_forces,face_weights);
			swimmer.particles.F(i)=-f;}
		iterate_face(axis,iter,fluid.mac_grid){const VectorDi& face=iter.Coord();
			if(face_weights(axis,face)!=(real)0)coupling_forces(axis,face)/=face_weights(axis,face);}
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);

		{std::string file_name=frame_dir+"/particles";
		swimmer.particles.Write_To_File_3d(file_name);}
		
		{std::string file_name=frame_dir+"/segment_mesh";
		swimmer.mesh->Write_To_File_3d(file_name);}

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

		//////Write fluid type
		//{std::string file_name=frame_dir+"/fluid_type";
		//Field<real,d> fluid_type;fluid_type.Resize(fluid.mac_grid.grid.cell_counts);
		//iterate_cell(iter,fluid.mac_grid.grid){const VectorDi& cell=iter.Coord();
		//	if(fluid.Is_Fluid_Cell(cell))fluid_type(cell)=(real)0;
		//	else fluid_type(cell)=(real)1;}
		//fluid_type.Write_To_File_3d(file_name);}

		//////Write interface
		//if(fluid.use_interface){
		//	std::string file_name=frame_dir+"/phi";
		//	fluid.levelset.phi.Write_To_File_3d(file_name);

		//	MarchingCubes<d> marching_cubes(fluid.levelset);
		//	marching_cubes.Marching();
		//	if(d==2){
		//		std::string file_name=frame_dir+"/segment_mesh";
		//		(*marching_cubes.mesh).Write_To_File_3d(file_name);}
		//	else{
		//		std::string file_name=frame_dir+"/triangle_mesh";
		//		(*marching_cubes.mesh).Write_To_File_3d(file_name);}}
	}
	
	virtual void Initialize()
	{
		swimmer.Initialize();

		int s=scale;real length=swimmer.L*((real)s/(real)swimmer.particles.Size());
		VectorDi cell_counts=VectorDi::Ones();cell_counts[0]=s;
		cell_counts[1]=cell_counts[0]/2;if(d>2)cell_counts[2]=cell_counts[0]/2;
		real dx=(real)length/cell_counts[0];
		VectorD domain_min=-dx*cell_counts.cast<real>()*(real).5;

		switch(test){
		case 1:{
			fluid.Initialize(cell_counts,dx,domain_min);
			//fluid.velocity.Fill((real).2,0);		
		}break;
		}

		coupling_forces.Resize(fluid.mac_grid.grid.cell_counts,(real)0);
	}

		//switch(test){
		//case 1:{
		//
		//}break;
		////case 1:{
		////	fluid.Initialize(cell_counts,(real)length/cell_counts[0]);
		////	fluid.velocity.Fill((real).2,0);

		////	////Source
		////	real source_speed=(real)1;int axis=0;
		////	iterate_face_in_one_dim(axis,iter,fluid.mac_grid){const VectorDi& face=iter.Coord();
		////		if(face[0]==0||face[0]==fluid.mac_grid.face_grids[0].node_counts[0]-1)fluid.Set_Psi_N(0,face,source_speed);}
		////	////Wall
		////	axis=1;iterate_face_in_one_dim(axis,iter,fluid.mac_grid){const VectorDi& face=iter.Coord();
		////		if(face[axis]==0||face[axis]==fluid.mac_grid.face_grids[axis].node_counts[axis]-1)fluid.Set_Psi_N(axis,face,0);}
		////	////Solid
		////	VectorD box_center=(real).5*(fluid.mac_grid.grid.domain_min+fluid.mac_grid.grid.domain_max);
		////	VectorD box_size=VectorD::Ones()*(real).1;Box<d> box(box_center-box_size,box_center+box_size);
		////	iterate_cell(iter,fluid.mac_grid.grid){const VectorDi& cell=iter.Coord();const VectorD& pos=fluid.mac_grid.grid.Center(cell);
		////		if(box.Inside(pos)){fluid.Set_Psi_D(cell,FluidEuler<d>::CellType::Solid);
		////			for(int axis=0;axis<d;axis++)for(int side=0;side<2;side++){VectorDi face=fluid.mac_grid.Cell_Incident_Face(axis,cell,side);fluid.Set_Psi_N(axis,face,(real)0);}}}
		////}break;
		////case 2:{
		////	fluid.Initialize(cell_counts,(real)length/cell_counts[0]);
		////	fluid.velocity.Fill((real).2,0);

		////	real source_speed=(real)1;
		////	VectorD box_center=fluid.mac_grid.grid.Center();
		////	box_center[0]=(real).2;VectorD box_size=VectorD::Ones()*(real).1;box_size[0]=fluid.mac_grid.grid.dx*(real)2;
		////	Box<d> box(box_center-box_size,box_center+box_size);

		////	////Initialize source
		////	VectorD source_velocity=VectorD::Unit(0)*(real)1;
		////	iterate_cell(iter,fluid.mac_grid.grid){const VectorDi& cell=iter.Coord();
		////		const VectorD& pos=fluid.mac_grid.grid.Center(cell);
		////			if(box.Inside(pos)){fluid.Set_Psi_D(cell,FluidEuler<d>::CellType::Solid);
		////			for(int axis=0;axis<d;axis++)for(int side=0;side<2;side++){
		////				VectorDi face=fluid.mac_grid.Cell_Incident_Face(axis,cell,side);fluid.Set_Psi_N(axis,face,source_velocity[axis]);}}}
		////}break;
		//}
		//fluid.Enforce_Boundary_Conditions();}
		// 
};

#endif