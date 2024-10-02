//#####################################################################
// Optimizer driver
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __OptimizerDriver_h__
#define __OptimizerDriver_h__
#include "Common.h"
#include "OptimizerIpOpt.h"
#include "OptimizerNlOpt.h"
#include "OptimizerSnOpt.h"
#include "OptimizerMma.h"

#ifdef USE_IPOPT
class SimpleOptimizerIpOpt : public OptimizerIpOpt
{
public:
	virtual void Initialize_Optimizer()
	{
		n_var=1024*64;
		n_hess=/*n_var*/0;
		Allocate_Data();

		for(int i=0;i<n_var;i++){var[i]=(real).0;var_L[i]=(real)-1.0;var_U[i]=(real)1000.;}
	}

	virtual real Compute_Objective(const real* var)
	{
		real f=(real)0;for(int i=0;i<n_var;i++)f+=pow(var[i]-(real)i,2);/*std::cout<<"f: "<<f<<std::endl;*/
		return f;
	}

	virtual void Compute_Gradient(const real* var,real* grad)
	{
		for(int i=0;i<n_var;i++){grad[i]=(real)2*(var[i]-(real)i);}
		//std::cout<<"Grad:\n";for(int i=0;i<n_var;i++)std::cout<<grad[i]<<", ";std::cout<<std::endl;
	}

	//virtual void Allocate_Hessian(int* iRow,int* jCol)
	//{
	//	for(int i=0;i<n_hess;i++){iRow[i]=i;jCol[i]=i;}
	//}

	//virtual void Compute_Hessian(real* hess_val)
	//{
	//	for(int i=0;i<n_hess;i++){hess_val[i]=(real)2;}
	//}

	virtual void Write_Substep(const int frame)
	{
		std::cout<<"---------- Frame "<<frame<<" ----------";std::cout<<std::endl;if(!write_intmed){return;}
		//std::cout<<"var: ";for(int i=0;i<n_var;i++)std::cout<<intmed_var[i]<<", ";std::cout<<std::endl;
		//std::cout<<"obj: "<<intmed_obj<<std::endl;
	}
};
#endif

#ifdef USE_NLOPT
class SimpleOptimizerNlOpt : public OptimizerNlOpt
{
public:
	virtual void Initialize_Optimizer()
	{
		n_var=5;
		Allocate_Data();

		for(int i=0;i<n_var;i++){var[i]=(real).0;var_L[i]=(real)-1.0;var_U[i]=(real)10.;}
	}

	virtual real Compute_Objective(const real* var)
	{
		real f=(real)0;for(int i=0;i<n_var;i++)f+=pow(var[i]-(real)i,2);/*std::cout<<"f: "<<f<<std::endl;*/
		return f;
	}

	virtual void Compute_Gradient(const real* var,real* grad)
	{
		for(int i=0;i<n_var;i++){grad[i]=(real)2*(var[i]-(real)i);}
		//std::cout<<"Grad:\n";for(int i=0;i<n_var;i++)std::cout<<grad[i]<<", ";std::cout<<std::endl;
	}

	virtual void Write_Substep(const int frame)
	{
		Seperation();std::cout<<"Frame "<<frame;Seperation();std::cout<<std::endl;
		std::cout<<"var: ";for(int i=0;i<n_var;i++)std::cout<<intmed_var[i]<<", ";std::cout<<std::endl;
		std::cout<<"obj: "<<intmed_obj<<std::endl;
	}
};
#endif

//class SimpleOptimizerSnOpt : public OptimizerSnOpt
//{typedef OptimizerSnOpt Base;typedef long long int TLI;
//public:
//	virtual void Initialize_Optimizer()
//	{
//		n_var=4;n_cons=1;nnz_cons_jac=n_var*n_cons;
//		Allocate_Data();
//		Initialize_Data();
//		Resize_Constraint_Jacobian();
//	}
//
//	virtual void Initialize_Data()
//	{
//		Base::Initialize_Data();
//		for(int i=0;i<n_var;i++){var_L[i]=(real)-1.;var_U[i]=(real)10.;var_state[i]=0;}
//	}
//
//	virtual real Compute_Objective(const real* var)
//	{
//		real f=(real)0;for(int i=0;i<n_var;i++)f+=pow(var[i]-(real)i,2);
//		return f;
//	}
//
//	virtual void Compute_Gradient(const real* var,real* grad)
//	{
//		for(int i=0;i<n_var;i++){grad[i]=(real)2*(var[i]-(real)i);}
//	}
//
//	virtual void Compute_Constraint(const real* var,real* constraint)
//	{
//		real c=(real)0;for(int i=0;i<n_var;i++){c+=var[i];}c-=(real)((n_var-1)*n_var)*(real).5;constraint[0]=c;
//	}
//
//	virtual void Resize_Constraint_Jacobian(TLI* row_num,TLI* col_num)
//	{
//		for(int i=0;i<n_var;i++){row_num[i]=1;col_num[i]=i;}
//	}
//	virtual void Resize_Constraint_Jacobian(){Base::Resize_Constraint_Jacobian();}
//
//	virtual void Compute_Constraint_Jacobian(const real* var,real* constraint_jacobian)
//	{
//		for(int i=0;i<n_var;i++){constraint_jacobian[i]=(real)1;}
//	}
//
//	virtual void Write_Substep(const int frame)
//	{
//		//Seperation();std::cout<<"Frame "<<frame;Seperation();std::cout<<std::endl;
//		//std::cout<<"var: ";for(int i=0;i<n_var;i++)std::cout<<intmed_var[i]<<", ";std::cout<<std::endl;
//		//std::cout<<"obj: "<<intmed_obj<<std::endl;
//	}
//};

class SimpleOptimizerMMA : public OptimizerMMA
{
public:
	virtual void Initialize_Optimizer()
	{
		n_var=10;
		n_cons=2;

		movlim=10;
		var_lb=-200.;
		var_ub=200.;

		Allocate_Data();

		for (int i=0;i<n_var;i++){
			var[i]=0.5;
			intmed_var[i]=0.5;}

		verbose=true;
	}

	virtual real Compute_Objective(const real* var)
	{
		real f=(real)0;for(int i=0;i<n_var;i++)f+=pow(var[i]-(real)i/(real)n_var,2);
		if(verbose){std::cout<<"comp obj: "<<f<<std::endl;}
		return f;
	}

	virtual void Compute_Gradient(const real* var,real* grad)
	{
		for(int i=0;i<n_var;i++){grad[i]=(real)2*(var[i]-(real)1/(real)n_var);}
		if(verbose){std::cout<<"comp grad: ";for(int i=0;i<n_var;i++)std::cout<<grad[i]<<", ";std::cout<<std::endl;}
	}

	virtual void Compute_Constraint(const real* var,real* constraint)
	{
		constraint[0]=0.0;
        for (int i=0;i<n_var;i++) {
            constraint[0] += var[i];
        }
        constraint[0]=constraint[0]/((real)n_var) - 1.0;

		constraint[1]=0.0;
        for (int i=0;i<n_var;i++) {
            constraint[1] += var[i];
        }
        constraint[1]=-(constraint[1]/((real)n_var) - 1.0);

		if(verbose){std::cout<<"comp cons: ";for(int i=0;i<n_cons;i++)std::cout<<constraint[i]<<", ";std::cout<<std::endl;}
	}

	virtual void Compute_Constraint_Grad(const real* var,real* constraint_grad)
	{
        for (int i=0;i<n_var;i++) {
            constraint_grad[n_cons*i+0]=1.0/((real)n_var);
            constraint_grad[n_cons*i+1]=-1.0/((real)n_var);}

		if(verbose){std::cout<<"comp cons grad: ";for(int i=0;i<n_cons;i++)for(int j=0;j<n_var;j++)std::cout<<constraint_grad[i*n_var+j]<<", ";std::cout<<std::endl;}
	}

	virtual void Write_Substep(const int frame)
	{
		std::cout<<"---------- Frame "<<frame<<" intmed var start ----------";std::cout<<std::endl;if(!write_intmed){return;}
		std::cout<<"substep var: ";for(int i=0;i<n_var;i++)std::cout<<intmed_var[i]<<", ";std::cout<<std::endl;
		std::cout<<"substep obj: "<<intmed_obj<<std::endl;
		std::cout<<"substep cons: ";for(int i=0;i<n_cons;i++)std::cout<<intmed_cons[i]<<", ";std::cout<<std::endl;
		std::cout<<"substep cons grad: ";for(int i=0;i<n_cons;i++)for(int j=0;j<n_var;j++)std::cout<<intmed_cons_grad[i*n_var+j]<<", ";std::cout<<std::endl;
		std::cout<<"---------- Frame "<<frame<<" intmed var end ----------";std::cout<<std::endl;
	}
};

class OptimizerTestDriver
{
public:
	int test;
	std::string output_dir;
	//SimpleOptimizerIpOpt opt;
	//SimpleOptimizerNlOpt opt;
	//SimpleOptimizerSnOpt opt;
	SimpleOptimizerMMA opt;

	void Initialize()
	{
		opt.Initialize_Optimizer();
	}

	void Run()
	{
		opt.Optimize();
	}
};
#endif
