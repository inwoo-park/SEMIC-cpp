#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"

#include "SemicParameters.h"
#include "SurfaceEnergyBalance.h"
#include "SemicArray.h"

namespace py = pybind11;

PYBIND11_MODULE(libpysemic, m){
	m.doc() = "SEMIC module in python.";

	py::class_<SemicParameters>(m, "Param")
	  .def(py::init<>())
	  .def_readwrite("ceff",&SemicParameters::ceff)
	  .def_readwrite("albi",&SemicParameters::albi)
	  .def_readwrite("albl",&SemicParameters::albl)
	  .def_readwrite("alb_smax",&SemicParameters::alb_smax)
	  .def_readwrite("alb_smin",&SemicParameters::alb_smin)
	  .def_readwrite("hcrit",&SemicParameters::hcrit)
	  .def_readwrite("rcrit",&SemicParameters::rcrit)
	  .def_readwrite("amp",&SemicParameters::amp)
	  .def_readwrite("csh",&SemicParameters::csh)
	  .def_readwrite("clh",&SemicParameters::clh)
	  .def_readwrite("tmin",&SemicParameters::tmin)
	  .def_readwrite("tmax",&SemicParameters::tmax)
	  .def_readwrite("tstic",&SemicParameters::tstic)
	  .def_readwrite("tsticsub",&SemicParameters::tsticsub)
	  .def_readwrite("tau_a",&SemicParameters::tau_a)
	  .def_readwrite("tau_f",&SemicParameters::tau_f)
	  .def_readwrite("w_crit",&SemicParameters::w_crit)
	  .def_readwrite("mcrit",&SemicParameters::mcrit)
	  .def_readwrite("tmid",&SemicParameters::tmid);

	py::class_<SemicConstants>(m, "Const") /*{{{*/
	  .def(py::init<>())
	  .def_readwrite("cls",&SemicConstants::cls)
	  .def_readwrite("clm",&SemicConstants::clm)
	  .def_readwrite("clv",&SemicConstants::clv)
	  .def_readwrite("cap",&SemicConstants::cap)
	  .def_readwrite("rhow",&SemicConstants::rhow)
	  .def_readwrite("hsmax",&SemicConstants::hsmax)
	  .def("Display",&SemicConstants::Display);
	/*}}}*/

	py::class_<SemicForcings>(m, "SemicForcings")
		.def(py::init<int, int>())
		.def_readwrite("nx",&SemicForcings::nx)
		.def_readwrite("ntime",&SemicForcings::ntime)
		.def_readwrite("sf",&SemicForcings::sf)
		.def_readwrite("rf",&SemicForcings::rf)
        .def_readwrite("t2m",&SemicForcings::t2m)
		.def_readwrite("sp",&SemicForcings::sp)
		.def_readwrite("lwd",&SemicForcings::lwd)
		.def_readwrite("swd",&SemicForcings::swd)
		.def_readwrite("wind",&SemicForcings::wind)
        .def_readwrite("rhoa",&SemicForcings::rhoa)
        .def_readwrite("qq",&SemicForcings::qq);	

	py::class_<SEMIC>(m, "SEMIC")
	   .def(py::init<>())
		.def_readwrite("Param",&SEMIC::Param)
		.def_readwrite("Const",&SEMIC::Const)
		/* initialize forcing variables */
		.def_readwrite("sf",&SEMIC::sf)
		.def_readwrite("rf",&SEMIC::rf)
		.def_readwrite("sp",&SEMIC::sp)
		.def_readwrite("lwd",&SEMIC::lwd)
		.def_readwrite("swd",&SEMIC::swd)
		.def_readwrite("wind",&SEMIC::wind)
		.def_readwrite("rhoa",&SEMIC::rhoa)
		.def_readwrite("qq",&SEMIC::qq)
		.def_readwrite("t2m",&SEMIC::t2m)
		/* model information */
		.def_readwrite("mask",&SEMIC::mask)
		.def_readwrite("alb_scheme",&SEMIC::alb_scheme)
		/* model results */
		.def_readwrite("tsurf",&SEMIC::tsurf)
		.def_readwrite("hsnow",&SEMIC::hsnow)
		.def_readwrite("hice",&SEMIC::hice)
		.def_readwrite("alb",&SEMIC::alb)
		.def_readwrite("alb_snow",&SEMIC::alb_snow)
		.def_readwrite("smb",&SEMIC::smb)
		.def_readwrite("smb_snow",&SEMIC::smb_snow)
		.def_readwrite("acc",&SEMIC::acc)
		.def_readwrite("melt",&SEMIC::melt)
		.def_readwrite("lwup",&SEMIC::lwup)
		.def_readwrite("shf",&SEMIC::shf)
		.def_readwrite("lhf",&SEMIC::lhf)
		.def_readwrite("qmr",&SEMIC::qmr)
		.def_readwrite("qmr_res",&SEMIC::qmr_res)
		.def_readwrite("verbose",&SEMIC::verbose)
		.def_readwrite("num_threads",&SEMIC::num_threads)
		.def("Initialize",&SEMIC::Initialize)
		.def("Display",&SEMIC::Display)
		.def("Albedo_Slater",&SEMIC::Albedo_Slater)
		.def("Albedo_Denby",&SEMIC::Albedo_Denby)
		.def("Albedo_ISBA",&SEMIC::Albedo_ISBA)
		.def("LongwaveRadiationUp",&SEMIC::LongwaveRadiationUp)
		.def("RunEnergyBalance",&SEMIC::RunEnergyBalance)
		.def("RunMassBalance",&SEMIC::RunMassBalance)
		.def("RunEnergyAndMassBalance",&SEMIC::RunEnergyAndMassBalance)
		.def("SetOpenmpThreads", (void (SEMIC::*)()) &SEMIC::SetOpenmpThreads)
		.def("SetOpenmpThreads", (void (SEMIC::*)(int)) &SEMIC::SetOpenmpThreads)
		.def("GetOpenmpThreads", &SEMIC::GetOpenmpThreads);

	py::class_<DoubleMatrix>(m, "DoubleMatrix")
        .def(py::init<size_t, size_t>())
		.def(py::init<py::array_t<double>>())
        .def("rows", &DoubleMatrix::rows)
        .def("cols", &DoubleMatrix::cols)
        .def("__call__", static_cast<const double& (DoubleMatrix::*)(size_t, size_t) const>(&DoubleMatrix::operator()))
        .def("__call__", static_cast<double& (DoubleMatrix::*)(size_t, size_t)>(&DoubleMatrix::operator()))
        .def("__str__", &DoubleMatrix::toString)
        .def("__repr__", &DoubleMatrix::toString)
		.def("to_numpy", &DoubleMatrix::toNumpy)
        .def("__getitem__", [](const DoubleMatrix& m, std::pair<size_t, size_t> index) {
            return m(index.first, index.second);
        })
        .def("__setitem__", [](DoubleMatrix& m, std::pair<size_t, size_t> index, double value) {
            m(index.first, index.second) = value;
        })
		.def("__getitem__", [](const DoubleMatrix& m, py::tuple index) -> py::object {
            if (index.size() != 2) throw std::invalid_argument("Index must be a tuple of size 2");

            py::object row_obj = index[0], col_obj = index[1];
            if (py::isinstance<py::slice>(row_obj) && py::isinstance<py::slice>(col_obj)) {
                throw std::invalid_argument("Slicing both rows and columns is not supported yet");
            }
            else if (py::isinstance<py::slice>(row_obj)) {
                py::slice row_slice = row_obj.cast<py::slice>();
                py::ssize_t start, stop, step, slicelength;
                if (!row_slice.compute(m.rows(), &start, &stop, &step, &slicelength))
                    throw py::error_already_set();
                std::vector<double> result(slicelength);
                for (py::ssize_t i = 0; i < slicelength; ++i) {
                    result[i] = m(start + i * step, col_obj.cast<py::ssize_t>());
                }
                return py::cast(result);
            }
            else if (py::isinstance<py::slice>(col_obj)) {
                py::slice col_slice = col_obj.cast<py::slice>();
                py::ssize_t start, stop, step, slicelength;
                if (!col_slice.compute(m.cols(), &start, &stop, &step, &slicelength))
                    throw py::error_already_set();
                std::vector<double> result(slicelength);
                for (py::ssize_t i = 0; i < slicelength; ++i) {
                    result[i] = m(row_obj.cast<py::ssize_t>(), start + i * step);
                }
                return py::cast(result);
            }
            else {
                return py::cast(m(row_obj.cast<py::ssize_t>(), col_obj.cast<py::ssize_t>()));
            }
        })
		.def("__getitem__", [](const DoubleMatrix& m, py::slice slice) -> py::object {
            py::ssize_t start, stop, step, slicelength;
            if (!slice.compute(m.rows(), &start, &stop, &step, &slicelength))
                throw py::error_already_set();
            std::vector<std::vector<double>> result(slicelength, std::vector<double>(m.cols()));
            for (py::ssize_t i = 0; i < slicelength; ++i) {
                for (size_t j = 0; j < m.cols(); ++j) {
                    result[i][j] = m(start + i * step, j);
                }
            }
            return py::cast(result);
        })
		.def("__setitem__", [](DoubleMatrix& m, py::tuple index, py::array_t<double> value) {
            if (index.size() != 2)
                throw std::invalid_argument("Index must be a tuple of size 2");

            py::object row_obj = index[0], col_obj = index[1];
            if (py::isinstance<py::slice>(row_obj) && py::isinstance<py::slice>(col_obj)) {
                throw std::invalid_argument("Slicing both rows and columns is not supported yet");
            }
            else if (py::isinstance<py::slice>(row_obj)) {
                py::slice row_slice = row_obj.cast<py::slice>();
                py::ssize_t start, stop, step, slicelength;
                if (!row_slice.compute(m.rows(), &start, &stop, &step, &slicelength))
                    throw py::error_already_set();
                if (slicelength != value.size())
                    throw std::invalid_argument("Slice length and input array size must match");
                auto buf = value.request();
                double* ptr = static_cast<double*>(buf.ptr);
                for (py::ssize_t i = 0; i < slicelength; ++i) {
                    m(start + i * step, col_obj.cast<py::ssize_t>()) = ptr[i];
                }
            }
            else if (py::isinstance<py::slice>(col_obj)) {
                py::slice col_slice = col_obj.cast<py::slice>();
                py::ssize_t start, stop, step, slicelength;
                if (!col_slice.compute(m.cols(), &start, &stop, &step, &slicelength))
                    throw py::error_already_set();
                if (slicelength != value.size())
                    throw std::invalid_argument("Slice length and input array size must match");
                auto buf = value.request();
                double* ptr = static_cast<double*>(buf.ptr);
                for (py::ssize_t i = 0; i < slicelength; ++i) {
                    m(row_obj.cast<py::ssize_t>(), start + i * step) = ptr[i];
                }
            } else {
                throw std::invalid_argument("Only row or column slice assignment is supported");
            }
        })
        .def("__setitem__", [](DoubleMatrix& m, py::slice slice, py::array_t<double> value) {
            py::ssize_t start, stop, step, slicelength;
            if (!slice.compute(m.rows(), &start, &stop, &step, &slicelength))
                throw py::error_already_set();
            if (slicelength != value.shape(0) || value.shape(1) != m.cols())
                throw std::invalid_argument("Slice length and input array dimensions must match");
            auto buf = value.request();
            double* ptr = static_cast<double*>(buf.ptr);
            for (py::ssize_t i = 0; i < slicelength; ++i) {
                for (size_t j = 0; j < m.cols(); ++j) {
                    m(start + i * step, j) = ptr[i * m.cols() + j];
                }
            }
        })
        .def("__setitem__", [](DoubleMatrix& m, py::slice slice, py::object value) {
            if (py::isinstance<py::array_t<double>>(value)) {
                py::array_t<double> array = value.cast<py::array_t<double>>();
                m.setMatrix(array);
            } else {
                throw std::invalid_argument("Only numpy arrays can be assigned to slices");
            }
        });
}
