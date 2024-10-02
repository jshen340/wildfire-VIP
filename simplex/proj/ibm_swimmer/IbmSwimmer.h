#ifndef __SoftSwimmer_h__
#define __SoftSwimmer_h__
#include "SoftBodyMassSpring.h"
#include "Particles.h"
#include "Constants.h"
#include "Mesh.h"

template<int d> class SoftSwimmer
{Typedef_VectorDii(d);
public:
	real L=(real)1;
	real a0=(real).01*L;
	real a1=(real).6*L;
	real P=(real)(.1*two_pi);
	real phase=(real)0;
	real mass=(real)100;

	VectorD mass_center=VectorD::Zero();
	VectorD linear_velocity=VectorD::Zero();
	Array<VectorD> objX;		////object space X
	Particles<d> particles;		////world space X
	std::shared_ptr<SegmentMesh<d> > mesh=nullptr;

	void Initialize()
	{
		int n=8;real step=L/(real)n;
		particles.Resize(n);
		for(int i=0;i<particles.Size();i++){
			particles.X(i)=step*(real)i*VectorD::Unit(0);
			particles.M(i)=mass/(real)n;
			particles.V(i)=VectorD::Zero();}
		mesh.reset(new SegmentMesh<d>(particles.XPtr()));
		for(int i=0;i<n-1;i++){mesh->elements.push_back(Vector2i(i,i+1));}

		mass_center=VectorD::Zero();
		for(int i=0;i<particles.Size();i++){
			mass_center+=particles.X(i);}
		mass_center/=(real)n;

		objX.resize(n);
		for(int i=0;i<n;i++){objX[i]=particles.X(i)-mass_center;}
	}

	void Advance(const real dt,const real time)
	{
		////mass center
		VectorD net_force=VectorD::Zero();
		for(int i=0;i<particles.Size();i++){
			net_force+=particles.F(i);}
		std::cout<<"net force: "<<net_force.transpose()<<std::endl;
		linear_velocity+=net_force/mass*dt;
		mass_center+=linear_velocity*dt;

		////obj space X
		for(int i=0;i<particles.Size();i++){
			real x=objX[i][0];
			real y=A(x)*sin(two_pi*(time/P-x/L)+phase);
			objX[i][1]=y;}
		
		////world space X
		for(int i=0;i<particles.Size();i++){
			VectorD old_pos=particles.X(i);
			particles.X(i)=mass_center+objX[i];
			particles.V(i)=(particles.X(i)-old_pos)/dt;}
	}

protected:
	real A(const real x) const {return a0+(a1-a0)*(x/L);}
};

#endif