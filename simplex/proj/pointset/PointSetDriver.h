#ifndef __PointSetDriver_h__
#define __PointSetDriver_h__

#include "PointSetDriver.h"
#include "Driver.h"
#include "PointSet.h"
#include "PointSetFunc.h"
#include "Mesh.h"
#include "AnalyticalFields.h"
#include "LeastSquares.h"
#include "RandomNumber.h"

using namespace AuxFunc;

template<int d> class PointSetDriver : public Driver
{Typedef_VectorDii(d);Typedef_MatrixD(d);using Base=Driver;
public:
	PointSet<d> ps;
	PointSet<d> sample_ps;
	VorticityField<d> field;
	GeometryParticles<d> aux_points;

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		switch(test){
		case 1:
		case 2:
		case 4:{
			for(int i=0;i<ps.points->Size();i++){
				ps.points->V(i)=field.Velocity(ps.points->X(i),time);
				ps.points->X(i)+=ps.points->V(i)*dt;}
			ps.Update();
			//ps.Update_Local_Frame(dt);
			ps.Reinitialize_Local_Frames();	
			ps.Point_Reseeding();
			//ps.Point_Relaxation();
		}break;
		case 3:{
			for(int i=0;i<ps.points->Size();i++){
				ps.points->V(i)=field.Velocity(ps.points->X(i),time);
				ps.points->X(i)+=ps.points->V(i)*dt;}
			ps.Update();
			//ps.Update_Local_Frame(dt);
			ps.Reinitialize_Local_Frames();		

			for(int i=0;i<aux_points.Size();i++){
				aux_points.X(i)=ps.Project_To_Surface(aux_points.X(i));}
		}break;
		case 5:{
			real rescale=.1;
			for(int i=0;i<ps.points->Size();i++){
				ps.points->X(i)[1]=rescale*sin(ps.points->X(i)[0]/rescale+time*(real)1);}
			ps.Update();
			ps.Reinitialize_Local_Frames();	
			ps.Point_Reseeding();
		}break;
		case 6: {
		}break;
		}
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);

		{std::string file_name=frame_dir+"/tracker_circles";
		PointSetFunc::Write_Tracker_Circles_To_File<d>(file_name,*ps.points);}

		{std::string file_name=frame_dir+"/segment_mesh";
		PointSetFunc::Write_Local_Frames_To_File(file_name,*ps.points);}
		
		if(aux_points.Size()!=0)
		{std::string file_name=frame_dir+"/tracker_points";
		aux_points.Write_To_File_3d_Fast(file_name);}

		{std::string file_name = frame_dir + "/point_force";
		Write_Segments_To_File_3d_Fast<d, real>(ps.points->XRef(), ps.points->FRef(), file_name); }
		{std::string file_name = frame_dir + "/point_velocity";
		Write_Segments_To_File_3d_Fast<d, real>(ps.points->XRef(), ps.points->VRef(), file_name); }

		std::cout<<"Write to frame "<<frame<<std::endl;
	}

	virtual void Initialize()
	{
		frame_rate=25;

		switch(test){
		case 1:{	////2D circle and 3D sphere, rigid rotation
			if constexpr (d==2){
				real dx=PointSetFunc::Initialize_Circle_Points(VectorD::Zero(),(real)1,64,*ps.points);
				ps.Initialize(dx,2);}
			else if constexpr (d==3){
				real dx=PointSetFunc::Initialize_Sphere_Points(VectorD::Zero(),(real)1,3,*ps.points);
				ps.Initialize(dx,2);}
			field.test=1;
		}break;
		case 2:{	////vortex motion
			VectorD c=VectorD::Unit(1)*(real).5;real R=(real).2;
			if constexpr (d==2){
				real dx=PointSetFunc::Initialize_Circle_Points(c,R,64,*ps.points);
				ps.Initialize(dx,2);}
			else if constexpr (d==3){
				real dx=PointSetFunc::Initialize_Sphere_Points(VectorD::Zero(),(real)1,3,*ps.points);
				ps.Initialize(dx,2);}
			field.test=2;
		}break;
		case 3:{	////test projection
			real dx;
			if constexpr (d==2){
				dx=PointSetFunc::Initialize_Circle_Points(VectorD::Zero(),(real)1,64,*ps.points);
				ps.Initialize(dx,2);}
			else if constexpr (d==3){
				dx=PointSetFunc::Initialize_Sphere_Points(VectorD::Zero(),(real)1,3,*ps.points);
				ps.Initialize(dx,2);}

			RandomNumber random((real)-1,(real)1);
			for(int i=0;i<ps.points->Size();i++){
				for(int k=0;k<4;k++){
					VectorD perturb=(real)2*dx*random.VectorValue<d>()+ps.points->X(i);
					int s=aux_points.Add_Element();
					aux_points.X(s)=perturb;}}
			field.test=2;
		}break;
		case 4:{	////test relaxation
			real dx;
			if constexpr (d==2){
				dx=PointSetFunc::Initialize_Circle_Points(VectorD::Zero(),(real)1,64,*ps.points);
				ps.Initialize(dx,2);}
			else if constexpr (d==3){
				dx=PointSetFunc::Initialize_Sphere_Points(VectorD::Zero(),(real)1,3,*ps.points);
				ps.Initialize(dx,2);}	
			ps.Update();

			////perturb the initial positions
			RandomNumber random((real)-1,(real)1);
			Array<VectorD> perturbed_x(ps.points->Size());
			for(int i=0;i<ps.points->Size();i++){
				Vector<real,d-1> t=(real)1.*dx*random.VectorValue<d-1>();
				VectorD p;ps.Unproject_To_World(t,ps.points->E(i),p);
				perturbed_x[i]=ps.points->X(i)+p;
				perturbed_x[i]=ps.Project_To_Surface(perturbed_x[i]);}
			for(int i=0;i<ps.points->Size();i++){
				ps.points->X(i)=perturbed_x[i];}
		}break;
		case 5:{	////sine wave
			real rescale=(real).1;
			int n=16;real length=two_pi*rescale;real dx=length/(real)n;
			if constexpr (d==2){
				PointSetFunc::Initialize_Segment_Points(VectorD::Zero(),VectorD::Unit(0)*length,n,*ps.points);
				ps.Initialize(dx,2);}
			else if constexpr (d==3){
				//PointSetFunc::Initialize_Rectangle_Points(n,n*2,dx,VectorD::Zero(),VectorD::Unit(1),*ps.points);
				PointSetFunc::Initialize_Circle_Points(VectorD::Unit(0)*length*(real).5,length*(real).5,VectorD::Unit(1),dx,*ps.points);
				ps.Initialize(dx,2);}
			ps.Update();

			std::cout<<"points: "<<ps.points->Size()<<std::endl;
			for(int i=0;i<ps.points->Size();i++){
				real y=rescale*sin(ps.points->X(i)[0]/rescale);
				ps.points->X(i)[1]=y;}
		}break;
		case 6:{	////test differential operator
			real dx;
			ps.dif_type=PointSet<d>::DifType::MLS;

			////init a unit sphere
			if constexpr (d==2){
				dx=PointSetFunc::Initialize_Circle_Points(VectorD::Zero(),(real)1,64,*ps.points);
				ps.Initialize(dx,2);}
			else if constexpr (d==3){
				dx=PointSetFunc::Initialize_Sphere_Points(VectorD::Zero(),(real)1,4,*ps.points);
				std::cout<<"Sample particles: dx="<<dx<<", pnum="<<ps.points->Size()<<std::endl;
				ps.Initialize(dx,2);}
			ps.Update();

			////sample particles on the sphere randomly
			if constexpr (d==2){////TOIMPL
			}
			else if constexpr (d==3){
				RandomNumber random(-half_pi,half_pi);
				Array<VectorD> rand_x(32);
				for(int i=0;i<rand_x.size();i++){
					real theta=random.Value()*2;
					real phi=random.Value();
					rand_x[i]=VectorD(sin(phi),cos(phi)*sin(theta),cos(phi)*cos(theta));
				}
				sample_ps.points->Add_Points(rand_x);
			}

			ps.Reinitialize_Local_Frames();
			ps.Update_Metric_Tensor();
			for(int i=0;i<ps.points->Size();i++){
				ps.points->M(i)=ps.points->X(i)[1];}

			////setup V=grad(m)
			std::cout<<"Compute grad(m) as points->V..."<<std::endl;
			std::function<real(const int)> Point_Mass=[&](const int idx)->real
				{return ps.points->M(idx);};
			for(int i=0;i<ps.points->Size();i++){
				VectorD grad_d=ps.Grad(i,Point_Mass);
				// grad_d=ps.points->E(i)*grad_d;
				ps.points->V(i)=grad_d;}
			
			////compute div(n)
			{
				real div_sum=0;
				std::cout<<"Compute div(n)..."<<std::endl;
				std::function<VectorD(const int)> Point_Normal=[&](const int idx)->VectorD
					{return ps.points->X(idx).normalized();};
				for(int i=0;i<ps.points->Size();i++){
					real div=ps.Div(i,Point_Normal);
					div_sum+=div;
				}
				std::cout<<"- avg curvature=div(n)="<<div_sum/ps.points->Size()<<std::endl;
			}

			{
				real div_sum=0;
				std::function<VectorD(const int)> Vector_One=[&](const int idx)->VectorD
					{return VectorD::Unit(0);};
				for(int i=0;i<ps.points->Size();i++){
					real div=ps.Div(i,Vector_One);
					div_sum+=div;
				}
				std::cout<<"- zero_div test: div(one)="<<div_sum/ps.points->Size()<<std::endl;
			}
			
			{
				////setup F=lap(x)/2
				real lap_sum=0;
				std::cout<<"Compute lap(x) as points->F..."<<std::endl;
				std::function<real(const int)> Point_Position_X=[&](const int idx)->real
					{return ps.points->X(idx)[0];};
				std::function<real(const int)> Point_Position_Y=[&](const int idx)->real
					{return ps.points->X(idx)[1];};
				std::function<real(const int)> Point_Position_Z=[&](const int idx)->real
					{return ps.points->X(idx)[2];};
				for(int i=0;i<ps.points->Size();i++){
					ps.points->F(i)[0]=ps.Laplacian(i,Point_Position_X);
					ps.points->F(i)[1]=ps.Laplacian(i,Point_Position_Y);
					if constexpr (d==3) ps.points->F(i)[2]=ps.Laplacian(i,Point_Position_Z);
					lap_sum+=ps.points->F(i).norm();
					ps.points->F(i)=-0.5*ps.points->F(i);}
				std::cout<<"- avg curvature=|lap(x)|="<<lap_sum/ps.points->Size()<<std::endl;
			}

			{
				////setup F=lap(x)/2 using CurvatureVector
				real lap_sum=0;
				for(int i=0;i<ps.points->Size();i++){
					ps.points->F(i)=ps.Curvature_Vector(i);
					lap_sum+=ps.points->F(i).norm();
					ps.points->F(i)=-0.5*ps.points->F(i);}
				std::cout<<"- avg curvature=|curv_vec(x)|="<<lap_sum/ps.points->Size()<<std::endl;
			}

			{
				real lap_sum=0;
				std::function<real(const int)> Real_One=[&](const int idx)->real
					{return 1;};
				for(int i=0;i<ps.points->Size();i++){
					lap_sum+=abs(ps.Laplacian_MLS(i,Real_One));}
				std::cout<<"- zero_lap test: |lap(1)|="<<lap_sum/ps.points->Size()<<std::endl;
			}

			////Operators on randomly sampled points
			{
				std::cout<<"Compute grad(m) as sample_ps->V, lap(x) as sample_ps->F..."<<std::endl;
				real lap_sum=0;
				for(int i=0;i<sample_ps.points->Size();i++){
					const VectorD& pos=sample_ps.points->X(i);
					MatrixD lf=ps.Local_Frame(pos);

					////Gradient
					std::function<real(const int)> Point_Mass=[&](const int idx)->real
					{return ps.points->M(idx);};
					VectorD grad_d=ps.Grad(pos,lf,Point_Mass);
					sample_ps.points->V(i)=grad_d;

					////Laplacian
					////setup F=lap(x)/2 using CurvatureVector
					sample_ps.points->F(i)=ps.Curvature_Vector(pos,lf);
					lap_sum+=ps.points->F(i).norm();
					sample_ps.points->F(i)=-0.5*sample_ps.points->F(i);
					// std::cout<<"---- pos["<<pos.transpose()<<"], grad["<<grad_d.transpose()<<"], lap["<<sample_ps.points->F(i).transpose()<<"]"<<std::endl;
					VectorD coef=ps.Fit_Local_Geometry_MLS(pos, lf);
				}
				std::cout<<"- avg curvature=|curv_vec(x)|="<<lap_sum/sample_ps.points->Size()<<std::endl;
			}

			////Add samplke_ps to ps for visualization
			int ps_size=ps.points->Size();
			ps.points->Add_Points(sample_ps.points->XRef());
			for(int i=0;i<sample_ps.points->Size();i++){
				ps.points->M(i+ps_size)=sample_ps.points->M(i);
				ps.points->V(i+ps_size)=sample_ps.points->V(i);
				ps.points->F(i+ps_size)=sample_ps.points->F(i);
			}

			// for(int i=0;i<ps.points->Size();i++){
			// 	std::cout<<"["<<i<<"]"<<ps.t_nbs[i][0]<<" pos("<<ps.points->X(i).transpose()<<") grad("<<ps.points->V(i).transpose()<<") lplc("<<ps.points->F(i).transpose()<<")"<<ps.points->F(i).norm()<<"("<<ps.points->E(i).col(2).transpose()<<")"<<ps.points->G(i)<<std::endl;
			// }
		}break;
		}
	}
};

#endif