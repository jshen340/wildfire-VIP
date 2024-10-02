//////////////////////////////////////////////////////////////////////////
// Use sparse solver to solve Poisson systems
// The reason is that some Poisson systems has very few unknowns against a large computational field
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "SparsePoissonSolverCPX.h"


template class SparsePoissonSolverCPX<Scalar, 2>;
template class SparsePoissonSolverCPX<Scalar, 3>;