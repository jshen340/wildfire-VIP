//////////////////////////////////////////////////////////////////////////
// Advection driver
//////////////////////////////////////////////////////////////////////////
#ifndef __AdvectionDriver_h__
#define __AdvectionDriver_h__
#include "Driver.h"
#include "MacGrid.h"
#include "Field.h"
#include "FaceField.h"
#include "FaceField.h"
#include "AnalyticalFields.h"
#include "Advection.h"
#include "AuxFunc.h"

template<int d> class AdvectionDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;
public:
	using Base::scale;

	real length_scale=(real).025;
	std::shared_ptr<VorticityField<d> > velocity_field=nullptr;

	MacGrid<d> mac_grid;
	Field<real,d> value;
	FaceField<real,d> velocity;

	virtual void Initialize()
	{
		//cfl=(real)1e-4;
		velocity_field=std::make_shared<VorticityField<d> >();
		velocity_field->test=4;

		VectorDi cell_counts=VectorDi::Ones()*scale;
		VectorD domain_size=VectorD::Ones()*(real)2; real dx=domain_size[0]/(real)cell_counts[0];
		mac_grid.Initialize(cell_counts,dx,VectorD::Zero());

		value.Resize(cell_counts,(real)0);
		velocity.Resize(cell_counts);

		iterate_face(axis,iter,mac_grid){const VectorDi& face=iter.Coord();
			const VectorD pos=mac_grid.Face_Center(axis,face);
			const VectorD v=velocity_field->Velocity(pos);
			velocity(axis,face)=v[axis];}

		switch(test){
		case 1:{
			VectorD center=mac_grid.grid.Position(AuxFunc::V<d>(.5,.75));
			real R=(real).1;real coef=(real)1e3;
			iterate_cell(iter,mac_grid.grid){const VectorDi& cell=iter.Coord();
				real r=(mac_grid.grid.Center(cell)-center).norm();
				if(r<R)value(cell)=coef*pow(R-r,3);}
		}break;}
	}

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		Field<real,d> ghost_value=value;
		Advection::Semi_Lagrangian<real,d>(dt,velocity,mac_grid,ghost_value,mac_grid,value);
	}

	virtual void Write_Output_Files(const int frame)
	{
		Base::Write_Output_Files(frame);

		{std::string file_name=output_dir+"/0/last_frame.txt";
		File::Write_Text_To_File(file_name,std::to_string(frame));}

		if(frame==0){
			std::string file_name=frame_dir+"/grid";
			mac_grid.grid.Write_To_File_3d(file_name);
			std::cout<<"Write to file "<<file_name<<std::endl;}

		{std::string file_name=frame_dir+"/velocity";
		velocity.Write_To_File_3d(file_name);}

		{std::string file_name=frame_dir+"/phi";
		value.Write_To_File_3d(file_name);}	
	}
};

#endif