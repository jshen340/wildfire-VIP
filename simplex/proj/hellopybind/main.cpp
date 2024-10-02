#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py=pybind11;

int Add(int a,int b){return a+b;}

py::array_t<double> Transfer_Array_From_Cpp_To_Python(double a, double b) 
{
	auto v = new std::vector<double>({a+1.,b+1.});
	auto capsule = py::capsule(v, [](void *v) { 
		std::cout<<"delete vector "<< ((*static_cast<std::vector<double>*>(v))[0])<<std::endl;delete reinterpret_cast<std::vector<int>*>(v); });
	return py::array_t<double>(v->size(), v->data(), capsule);
}

PYBIND11_MODULE(/*package name in Python*/hellopybind,m)
{
	m.def(/*function name in Python*/"add",/*function name in C++*/&Add,/*comments*/"hellopybind add");
	m.def(/*function name in Python*/"transfer",/*function name in C++*/&Transfer_Array_From_Cpp_To_Python,/*comments*/"hellopybind transfer");
}