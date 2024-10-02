#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "Common.h"

namespace py = pybind11;

int main() 
{
	////pip install numpy scipy
	py::scoped_interpreter guard{}; // start the interpreter and keep it alive
	py::print("Hello, pybind!"); // use the Python API

	py::module np = py::module::import("numpy");
	py::module scipy = py::module::import("scipy.linalg");
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3, 3);
	A << 0, 0.5, 0.5, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4;
	std::cout << "A: \n" << A << std::endl;
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 2);
	B << 0.3, 0.4, 0.1, 0.2, 0.2, 0.2;
	std::cout << "B: \n" << B << std::endl;
	Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(3, 3);
	std::cout << "Q: \n" << Q << std::endl;
	Eigen::MatrixXd R = Eigen::MatrixXd::Identity(2, 2);
	std::cout << "R: \n" << R << std::endl;

	py::function solve = scipy.attr("solve_continuous_are");
	py::object retVals = solve(A, B, Q, R);

	Eigen::MatrixXd result = retVals.cast<Eigen::MatrixXd>();
	std::cout << "printing in c++\n" << result << std::endl;
	std::cout<<"successfully finish pybind calls"<<std::endl;
}