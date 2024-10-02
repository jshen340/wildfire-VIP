#pragma once
#include "para.h"

class LinearMapping
{
public:
	virtual int xDoF() = 0;

	virtual int yDoF() = 0;

	//input p, get Ap
	virtual void applyMapping(Scalar *Ap, Scalar *p) = 0;

	virtual void update(){}
};