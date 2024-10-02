#include <iostream>
#include "AutoGradForwardOptimizer.h"

#ifndef __Main_cpp__
#define __Main_cpp__

int main(int argc, char* argv[])
{
	AutoGradForwardOptimizer test_optimizer;
	test_optimizer.Initialize_Optimizer();
	test_optimizer.Optimize();
}

#endif