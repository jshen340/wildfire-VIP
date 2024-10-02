#pragma once
#include "Common.h"
#include "File.h"
#include "GeometryPrimitives.h"
#include "FluidParticleLevelset.h"
#include "Particles.h"
#include "MarchingCubes.h"
#include "Driver.h"

template<int d> class ParticleLevelsetAdvectionDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;
public:
	FluidParticleLevelset<d> fluid;

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		
        Set_Zalesak_Velocity(fluid.mac_grid,fluid.velocity);
		// fluid.Advance(dt);
        fluid.Advection(dt);

		if(int(floor(time/dt))%20==0){
			fluid.pls.Reseed();//// TOFIX: remove frame counter
			fluid.pls.Update_Particle_Radii();
		}

        //// reset velocity field
        Set_Zalesak_Velocity(fluid.mac_grid,fluid.velocity);
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

		////Write interface
		std::string file_name=frame_dir+"/phi";
		fluid.levelset.phi.Write_To_File_3d(file_name);

		MarchingCubes<d> marching_cubes(fluid.levelset);
		marching_cubes.Marching();
		if(d==2){
			std::string file_name=frame_dir+"/segment_mesh";
			(*marching_cubes.mesh).Write_To_File_3d(file_name);}
		else{
			std::string file_name=frame_dir+"/triangle_mesh";
			(*marching_cubes.mesh).Write_To_File_3d(file_name);}

        {std::string file_name=frame_dir+"/particles";
		    fluid.pls.particles.Write_To_File_3d(file_name);}

		{std::string file_name = frame_dir + "/point_impulse";
		for(int i=0;i<fluid.pls.particles.Size();i++){
			VectorD norm=fluid.pls.Normal(fluid.pls.particles.X(i));
			real r=fluid.pls.Particle_Radius(i);
			fluid.pls.particles.F(i)=-fluid.pls.Particle_Sign(i)*norm.normalized()*r;}
		Write_Segments_To_File_3d_Fast<d, real>(fluid.pls.particles.XRef(), fluid.pls.particles.FRef(), file_name);
		}

		{std::string file_name = frame_dir + "/point_velocity";
		Write_Segments_To_File_3d_Fast<d, real>(fluid.pls.particles.XRef(), fluid.pls.particles.VRef(), file_name);
		}
		
	}
	
	virtual void Initialize()
	{
		int s=scale;real length=(real)1;VectorDi cell_counts=VectorDi::Ones();cell_counts[0]=s;
		cell_counts[1]=cell_counts[0]/2;if(d>2)cell_counts[2]=cell_counts[0]/2;

		switch(test){
		case 1:{//// particle levelset
			cell_counts=VectorDi::Ones()*100;
			fluid.Initialize(cell_counts,(real)length/cell_counts[0]);
			fluid.use_body_force=false;
			
			Init_Zalesak_Sphere(fluid.mac_grid, fluid.levelset);
            fluid.Setup_Particle_Level_Set();
			fluid.projection.use_multigrid_solver=true;
		}break;
		case 2:{//// standard levelset
			fluid.use_particle=false;
			cell_counts=VectorDi::Ones()*100;
			fluid.Initialize(cell_counts,(real)length/cell_counts[0]);
			fluid.use_body_force=false;
			
			Init_Zalesak_Sphere(fluid.mac_grid, fluid.levelset);
			fluid.projection.use_multigrid_solver=true;
		}break;}

		// fluid.Enforce_Boundary_Conditions();
		// fluid.Advance((real).01);
	}

    void Init_Zalesak_Sphere(const MacGrid<d>& mac_grid, LevelSet<d>& levelset) {
        if constexpr(d==2){
            const VectorD& length=mac_grid.grid.Length();
            const VectorD& domain_min=mac_grid.grid.domain_min;
            VectorD center=domain_min+length*0.5f+VectorD::Unit(1)*length(1)*0.25;
            real r=length(0)*0.15;
            VectorD bottom=center-r*VectorD::Unit(1);
            Sphere<d> sphere(center,r);
            //// slot: 5 width, 12.5 deep
            Box<d> box(bottom+VectorD(-0.025*length(0),0),bottom+VectorD(0.025*length(0),0.125*length(1)));
			std::cout<<"BOX:"<<box.min_corner.transpose()<<"]["<< box.max_corner.transpose()<<std::endl;
            iterate_cell(iter,fluid.mac_grid.grid){const VectorDi& cell=iter.Coord();
				VectorD pos=fluid.mac_grid.grid.Center(cell);
				// fluid.levelset.phi(cell)=box.Phi(pos);
                if(sphere.Phi(pos)>0)
                    fluid.levelset.phi(cell)=sphere.Phi(pos);
                else
					fluid.levelset.phi(cell)=sphere.Phi(pos);
                    fluid.levelset.phi(cell)=std::max(-box.Phi(pos),sphere.Phi(pos));
            }
        }
    }
	
    void Set_Zalesak_Velocity(const MacGrid<d>& mac_grid, FaceField<real, d>& velocity) {
		const VectorD& domain_min=mac_grid.grid.domain_min;
        iterate_face_in_one_dim(0,iter,mac_grid){
            VectorDi face = iter.Coord();
            VectorD face_center=mac_grid.Face_Center(0,face);
            velocity(0,face)=pi/314*(0.5*mac_grid.grid.Length()(1)-face_center(1)+domain_min(1));
        }

        iterate_face_in_one_dim(1,iter,mac_grid){
            VectorDi face = iter.Coord();
            VectorD face_center=mac_grid.Face_Center(1,face);
            velocity(1,face)=pi/314*(face_center(0)-domain_min(0)-0.5*mac_grid.grid.Length()(0));
        }

        if constexpr(d==3){
            iterate_face_in_one_dim(1,iter,mac_grid){
                VectorDi face = iter.Coord();
                velocity(2,face)=0;
            }
        }

    }
};