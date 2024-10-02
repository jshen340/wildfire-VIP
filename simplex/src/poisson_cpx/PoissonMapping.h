#pragma once
#include "LinearMapping.h"
#include "form01mapping.h"
#include "PoissonDescriptor.h"

class PoissonMapping :public LinearMapping
{
public:
	int Nx, Ny;
	int dof;
	grid2D grid;
	const PoissonDescriptor<2> *descr;

	Scalar *d_temp;

	void init(PoissonDescriptor<2> *_descr);

	void Destroy(void) {
		AuxFuncCPX::Global_Free(d_temp, DataHolder::DEVICE);
	}

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

class PoissonMappingFixed :public LinearMapping
{
public:
	int Nx, Ny;
	int dof;
	grid2D grid;
	const PoissonDescriptor<2> *descr;

	Scalar *d_temp;

	void init(PoissonDescriptor<2> *_descr);

	void Destroy(void) {
		AuxFuncCPX::Global_Free(d_temp, DataHolder::DEVICE);
	}

	int xDoF() override;

	int yDoF() override;
	//p->Ap
	void applyMapping(Scalar *Ap, Scalar *p) override;
};