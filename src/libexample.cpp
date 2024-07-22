// src/example.cpp
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "omp.h"
#include <iostream>
#include <vector>
#include "math.h"

namespace py = pybind11;
using namespace std;

int add(int i, int j) {
    return i + j;
}

/* Test pybind11 with openmp */
std::vector<double> multiply_by_two(const std::vector<double>& vec, int num_threads) {
	int nx = static_cast<int>(vec.size());
	std::vector<double> result(nx,0.0);

	omp_set_num_threads(num_threads);

	#pragma omp parallel
	{
		#pragma omp single
		{
			cout << "number of threads = " << omp_get_num_threads() << endl;
		}

		#pragma omp for
		for (int i=0; i<nx; i++){
			result[i] = sqrt(vec[i]);
		}
	}
	return result;
}

PYBIND11_MODULE(libexample, m) {
	  m.def("add", &add, "A function which adds two numbers").
	  def("multiply_by_two", &multiply_by_two, "Multiply each element of the vector by two.");
}
