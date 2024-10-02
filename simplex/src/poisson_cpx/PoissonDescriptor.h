#pragma once

#include "cpuUtils.h"
#include "form01mapping.h"
#include "grid3D.h"
#include "TypeFunc.h"
#include "AuxFuncCPX.h"

template<int d>
class PoissonDescriptor
{
	Typedef_VectorDii(d);
public:
	VectorDi N;
	//int Nx, Ny;
	int size; //dof
	int fsize;
	using Grid = typename If<d == 2, grid2D, grid3D >::Type;
	Grid grid;
	Scalar* h_vol = nullptr, * d_vol = nullptr;//h on CPU and d on GPU
	bool* h_fixed = nullptr, * d_fixed = nullptr;

public:


	void init(const VectorDi& _n);
	//void init(int _Nx, int _Ny);

	void Destroy(void) {
		//PoissonDescriptor is sometimes copied in Poisson
		//So we can't write a deconstructor
		AuxFuncCPX::Global_Free(h_vol, DataHolder::HOST);
		AuxFuncCPX::Global_Free(d_vol, DataHolder::DEVICE);
		AuxFuncCPX::Global_Free(h_fixed, DataHolder::HOST);
		AuxFuncCPX::Global_Free(d_fixed, DataHolder::DEVICE);
	}

	void setFixed(bool *fixed);

	void setVol(Scalar *vol);

	void finish();

	void toHost();

	void toDevice();

};