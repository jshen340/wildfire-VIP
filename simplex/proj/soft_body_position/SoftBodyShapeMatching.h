#ifndef __SoftBodyShapeMatching_h__
#define __SoftBodyShapeMatching_h__

#include "Particles.h"
#include "BoundaryCondition.h"
#include "AuxFunc.h"
using namespace AuxFunc;

template<int d> class SoftBodyShapeMatching
{Typedef_VectorDii(d);Typedef_MatrixD(d);
public:
	Particles<d> particles;
	Array<VectorD> q0;
	MatrixD A_qq=MatrixD::Zero();
	int n=0;
	real alpha=(real).9;
	real beta=(real).5;
	bool use_body_force=true;
	VectorD g=VectorD::Unit(1)*-1;
	std::function<void(const real dt)> Collision;
	std::function<void(const real dt,const real time)> Kinematic_Boundary_Condition;
	BoundaryConditionMesh<d> bc;

	void Initialize()
	{
		n=particles.Size();
		q0=particles.XRef();
		VectorD c0=Mean<VectorD>(q0);
		for(auto& q:q0)q=q-c0;
		A_qq=MatrixD::Zero();
		for(int i=0;i<n;i++){A_qq+=q0[i]*q0[i].transpose();}
		A_qq=A_qq.inverse().eval();
	}

	virtual void Advance(const real dt,const real time)
	{
		Fill(particles.FRef(),VectorD::Zero());
		if(use_body_force)for(int i=0;i<n;i++)particles.F(i)+=particles.M(i)*g;
		for(auto& iter:bc.forces){particles.F(iter.first)+=iter.second;}

		if(Collision!=nullptr)Collision(dt);
		if(Kinematic_Boundary_Condition!=nullptr)Kinematic_Boundary_Condition(dt,time);

		////shape matching
		VectorD c=Mean<VectorD>(particles.XRef());
		Array<VectorD> p0=particles.XRef();
		for(auto& p:p0)p=p-c;
		MatrixD A_pq=MatrixD::Zero();
		for(int i=0;i<n;i++){A_pq+=p0[i]*q0[i].transpose();}
		MatrixD R,S;Polar_Decomposition<d>(A_pq,R,S);
		MatrixD A=A_pq*A_qq;
		real det=A.determinant();
		A/=pow(det,one_third);
		MatrixD B=((real)1-beta)*A+beta*R;
		for(int i=0;i<n;i++){
			VectorD g=B*q0[i]+c;

			////add an extra damping force
			real damping=(real)2;
			VectorD xij=(g-particles.X(i)).normalized();
			VectorD f_damp=damping*(-particles.V(i)).dot(xij)*xij;
			particles.F(i)+=f_damp;

			particles.V(i)+=alpha*(g-particles.X(i))/dt;
			particles.V(i)+=particles.F(i)/particles.M(i)*dt;
			particles.X(i)+=particles.V(i)*dt;}
	}
};
#endif