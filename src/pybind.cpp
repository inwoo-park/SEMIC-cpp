#include <pybind11/pybind11.h>
#include "SurfaceEnergyBalance.h"

namespace py = pybind11;

PYBIND11_MODULE(semic, m){
	m.doc() = "SEMIC module in python.";

	py::class_<SEMIC>(m, "SEMIC")
		.def(py::init<>())
		.def("Initialize",&SEMIC::Initialize);
}