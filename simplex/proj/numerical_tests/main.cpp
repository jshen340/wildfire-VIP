#include <iostream>
#include "Integration.h"
#include "Kernels.h"
#include "RandomNumber.h"
#include "Timer.h"
//using namespace std;

void Test_Line_Int(void) {
	std::function<real(const Vector1&)> func = [&](const Vector1& x)->real {return pow(x[0], 9); };
	Vector1 a(0), b(1);
	std::cout << "Test integration in [0,1] with f=x^5\n Expected result: 0.1\n";
	for (int degree = 2; degree <= 4; degree++) {
		std::cout << "Numerical result with degree=" << degree << " and n=10: " << Integration::Newton_Cotes_Line<1, real>(a, b, func, 10, degree) << "\n";
		std::cout << "Numerical result with degree=" << degree << " and n=100: " << Integration::Newton_Cotes_Line<1, real>(a, b, func, 100, degree) << "\n";
	}
}

template<int d>
real Kernel_Integration(KernelType mode, KernelSPH kernel, int n) {
	Vector<real, d> min_crn = Vector<real, d>::Ones() * (-kernel.h);
	Vector<real, d> max_crn = Vector<real, d>::Ones() * (+kernel.h);
	Vector<int, d> res = Vector<int, d>::Ones() * n;
	std::function<real(const Vector<real, d>&)> f = [&](const Vector<real, d>& r)->real {
		return kernel.Weight(d, r.norm(), mode);
	};
	return Integration::Newton_Cotes<d, real>(min_crn, max_crn, f, res);
}

void Test_Kernels(int n) {
	real h = 5.0;
	KernelType modes[5] = { KernelType::POLY6, KernelType::SPIKY, KernelType::CUBIC, KernelType::QUINTIC, KernelType::GAUSSIAN };
	std::string names[5] = { "POLY6","SPIKY","CUBIC","QUINTIC","GAUSSIAN" };
	KernelSPH kernel(h);
	for (int i = 0; i < 5; i++) {
		std::cout << "Start to test kernel " << names[i] << "\n";
		std::cout << "1D integration: " << Kernel_Integration<1>(modes[i], kernel, n) << "\n";
		std::cout << "2D integration: " << Kernel_Integration<2>(modes[i], kernel, n) << "\n";
		std::cout << "3D integration: " << Kernel_Integration<3>(modes[i], kernel, n) << "\n";
		std::cout << "\n";
	}
}

template<int d>
void Point_Grad_Test(const Vector<real, d>& pos, const KernelSPH& kernel, KernelType mode, int n) {
	//Gradient
	Vector<real, d> grad_ana = kernel.Grad<d>(pos, mode), grad_num;
	for (int axis = 0; axis < d; axis++) {
		Vector<real, d> pos1 = pos, pos2 = pos;
		pos1[axis] -= kernel.h / n;
		pos2[axis] += kernel.h / n;
		real f1 = kernel.Weight<d>(pos1, mode), f2 = kernel.Weight<d>(pos2, mode);
		grad_num[axis] = (f2 - f1) / (kernel.h / n * 2);
	}
	Vector<real, d> grad_err = grad_num - grad_ana;
	std::cout << "Position: " << pos.transpose() << " ,error norm=" << grad_err.norm() << "\n";
}

//template<int d>
//void Point_Laplacian_Test(const Vector<real, d>& pos, const KernelSPH<d>& kernel, int n) {
//	KernelType mode = KernelType::POLY6;
//	real lap_ana = kernel.Lap_Poly6(pos), lap_num = 0;
//	real f0 = kernel.Weight(pos.norm(), mode);
//	for (int axis = 0; axis < d; axis++) {
//		Vector<real, d> pos1 = pos, pos2 = pos;
//		pos1[axis] -= kernel.h / n;
//		pos2[axis] += kernel.h / n;
//		real f1 = kernel.Weight(pos1.norm(), mode), f2 = kernel.Weight(pos2.norm(), mode);
//		lap_num += (f1 + f2 - 2 * f0) / pow(kernel.h / n, 2);
//	}
//	real err = lap_num - lap_ana;
//	std::cout << "Predicted: " << lap_ana << " ,numerical: " << lap_num << " ,error=" << err << "\n";
//}
//
//template<int d>
//void Test_Laplacian(real h, int n) {
//	KernelSPH<d> kernel(h);
//	int tests = 10;
//	for (int iter = 0; iter < tests; iter++) {
//		Vector<real, d> pos = RandomFunc::Random_Vector_Spherical<d>(h);
//		Point_Laplacian_Test<d>(pos, kernel, n);
//	}
//}

template<int d>
void Test_Differential(real h, int n) {
	KernelType modes[5] = { KernelType::POLY6, KernelType::SPIKY, KernelType::CUBIC, KernelType::QUINTIC, KernelType::GAUSSIAN };
	std::string names[5] = { "POLY6","SPIKY","CUBIC","QUINTIC","GAUSSIAN" };
	KernelSPH kernel(h);
	int tests = 20;
	for (int i = 0; i < 5; i++) {
		AuxFunc::Seperation_Line();
		std::cout << "Start to test kernel " << names[i] << "\n";
		for (int iter = 0; iter < tests; iter++) {
			Vector<real, d> pos = RandomFunc::Random_Vector_Spherical<d>(h);
			Point_Grad_Test<d>(pos, kernel, modes[i], n);
		}
		AuxFunc::Seperation_Line(); std::cout << "\n";
	}
}

void Performance_Test(void) {
	real A = 2.0;
	real B = 0.5;
	int N = 1e8;
	Vector3 vec_eg(9, 9, 8), vec_eg1;
	real a[3] = { 9,9,8 };
	Timer timer;
	for (int i = 0; i < N; i++) {
		Vector3 vec_temp(i, i, i);
		vec_eg += vec_temp;
	}
	std::cout << "vec_eg: " << vec_eg.transpose() << "\n";
	timer.Elapse_And_Output_And_Reset("Eigen vector operation");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < 3; j++) {
			a[j] *= A;
			a[j] *= B;
		}
	}
	std::cout << "a: " << a[0] << " " << a[1] << " " << a[2] << "\n";
	timer.Elapse_And_Output_And_Reset("Array operation");
}

int main() {
	//Performance_Test();
	//std::cout << "\n=================================\nTest Kernel Integrations:";
	Test_Kernels(100);

	real h = 1.0 / 16 * 3;
	std::cout << "h: " << h << "\n";
	std::cout << "\n=================================\nTest 1D Grad:";
	Test_Differential<1>(h, 1000);
	//std::cout << "\n=================================\nTest 2D Grad:";
	//Test_Differential<2>(h, 1000);
	//std::cout << "\n=================================\nTest 3D Grad:";
	//Test_Differential<3>(h, 1000);

	//Test_Laplacian<1>(5.0, 100);
	return 0;
}
