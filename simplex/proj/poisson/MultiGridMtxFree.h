//////////////////////////////////////////////////////////////////////////
// Multigrid for matrix-free Krylov systems
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __MultiGridMtxFree_h__
#define __MultiGridMtxFree_h__
#include "MultiGrid.h"
#include "Field.h"
#include "FaceField.h"
#include "BitField.h"
#include "BitFaceField.h"

namespace MultiGrid{

//////////////////////////////////////////////////////////////////////////
////Matrix-free multigrid data structure on CPU
template<int d,class T,class T_SPARSE,class T_VECTOR> class MultiGridSystemMtxFree
{public:
	Params params;
	Array<T_SPARSE*> A;
	Array<T_VECTOR*> b;
	Array<T_VECTOR*> x;
	Array<T_VECTOR*> r;
	Array<Vector<int,d> > coarsen_factor;

	virtual void V_Cycle()
	{
		//std::cout<<"x: ";
		//(*x[0]).Print();
		//std::cout<<"b: ";
		//(*b[0]).Print();

		//(*A[0]).Gauss_Seidel_Smoothing(*x[0],*b[0],40,false);	////

		////std::cout<<"x: ";
		////(*x[0]).Print();
		////std::cout<<"b: ";
		////(*b[0]).Print();

		for(int i=0;i<=params.levels-2;i++){
			////smooth A_i*x_i=b_i
			if(i>0)(*x[i]).Set_Zero();
			(*A[i]).Gauss_Seidel_Smoothing(*x[i],*b[i],params.smooth_num,false);
			////r_i=b_i-A_i*x_i
			(*A[i]).Res((*b[i]),(*x[i]),(*r[i]));
			////b_{i+1}=R_i*r_i
			(*r[i]).Restriction((*b[i+1]),coarsen_factor[i]);}

		////////solve A_{n-1}x_{n-1}=b_{n-1} with a direct solver
		//(*A[params.levels-1]).Solve((*x[params.levels-1]),(*b[params.levels-1]));
		int i=params.levels-1;
		(*A[i]).Gauss_Seidel_Smoothing(*x[i],*b[i],20,false);
		(*A[i]).Gauss_Seidel_Smoothing(*x[i],*b[i],20,true);

		for(int i=params.levels-2;i>=0;i--){
			////x_i=x_i+P_i*x_{i+1}
			(*x[i+1]).Prolongation(*r[i],coarsen_factor[i]);
			(*x[i]).Axpy((real)1,(*x[i]),(*r[i]));
			////smooth A_i*x_i=b_i
			(*A[i]).Gauss_Seidel_Smoothing(*x[i],*b[i],params.smooth_num,true);}

		////std::cout<<"x: ";
		////(*x[0]).Print();
	}
};

//////////////////////////////////////////////////////////////////////////
////These functions are performance-sensitive!

////restrict cells
template<int d,class T=real> void Restriction(const Field<T,d>& fine,Field<T,d>& coarse,const Vector<int,d>& factor)
{
	if constexpr (d==2){
		const T e1=(T)1/(T)8;
		const T e2=(T)3/(T)8;
		const T e3=(T)3/(T)8;
		const T e4=(T)1/(T)8;

		const T coef[4][4]={e1*e1,e1*e2,e1*e3,e1*e4,
							e2*e1,e2*e2,e2*e3,e2*e4,
							e3*e1,e3*e2,e3*e3,e3*e4,
							e4*e1,e4*e2,e4*e3,e4*e4};

		auto& coarse_count=coarse.counts;
		auto& fine_count=fine.counts;
		for(int i=0;i<coarse_count[0];i++){
			for(int j=0;j<coarse_count[1];j++){
				T c=(T)0;
				for(int ix=-1;ix<=2;ix++){
					for(int iy=-1;iy<=2;iy++){
						int x=i*factor[0]+ix;int y=j*factor[1]+iy;
						if(x<0||x>=fine_count[0]||y<0||y>=fine_count[1]){continue;}
						c+=fine(x,y)*coef[ix+1][iy+1];}}
				coarse(i,j)=c;}}	
	}
	else if constexpr (d==3){
		////TODO
	}
}

////prolongate cells
template<int d,class T=real> void Prolongation(const Field<T,d>& coarse,Field<T,d>& fine,const Vector<int,d>& factor)
{
	if constexpr (d==2){
		const T e1=(T)1/(T)8;
		const T e2=(T)3/(T)8;
		const T e3=(T)3/(T)8;
		const T e4=(T)1/(T)8;

		const T coef[4][4]={e1*e1,e1*e2,e1*e3,e1*e4,
							e2*e1,e2*e2,e2*e3,e2*e4,
							e3*e1,e3*e2,e3*e3,e3*e4,
							e4*e1,e4*e2,e4*e3,e4*e4};

		auto& coarse_count=coarse.counts;
		auto& fine_count=fine.counts;
		fine.Fill((T)0);
		real rescale=Pow((real)2,d);
		for(int i=0;i<coarse_count[0];i++){
			for(int j=0;j<coarse_count[1];j++){
				T c=coarse(i,j);
				for(int ix=-1;ix<=2;ix++){
					for(int iy=-1;iy<=2;iy++){
						int x=i*factor[0]+ix;int y=j*factor[1]+iy;
						if(x<0||x>=fine_count[0]||y<0||y>=fine_count[1]){continue;}
						fine(x,y)+=rescale*c*coef[ix+1][iy+1];}}}}	
	}
	else if constexpr (d==3){
		////TODO
	}
}

////restrict boundary conditions
template<int d,class T=real> void Restriction_Boundary_Condition(const BitField<d>& fine,BitField<d>& coarse,
	const Field<T,d>& fine_val,Field<T,d>& coarse_val,const Vector<int,d>& factor)
{
	if constexpr (d==2){
		auto& coarse_count=coarse.counts;
		auto& fine_count=fine.counts;
		for(int i=0;i<coarse_count[0];i++){
			for(int j=0;j<coarse_count[1];j++){
				T c=(T)0;int n=0;
				for(int ix=0;ix<factor[0];ix++){
					for(int iy=0;iy<factor[1];iy++){
						int x=i*factor[0]+ix;int y=j*factor[1]+iy;
						if(fine(x,y)){c+=fine_val(x,y);n++;}}}
				if(n>0){coarse_val(i,j)=c/(T)n;coarse.Set(i,j,true);}
				else {coarse_val(i,j)=(T)0;coarse.Set(i,j,false);}}}			
	}
	else if constexpr (d==3){
	
	}
}

template<int d,class T=real> void Restriction_Boundary_Condition(const BitFaceField<d>& fine,BitFaceField<d>& coarse,
	const FaceField<T,d>& fine_val,FaceField<T,d>& coarse_val,const Vector<int,d>& factor)
{
	////TODO		
}

};
#endif
