#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "SemicParameters.h"
#include "SurfaceEnergyBalance.h"

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
		.def("Initialize",&SEMIC::Initialize)
		.def("Display",&SEMIC::Display)
		.def("Albedo_Slater",&SEMIC::Albedo_Slater)
		.def("Albedo_Denby",&SEMIC::Albedo_Denby)
		.def("Albedo_ISBA",&SEMIC::Albedo_ISBA)
		.def("LongwaveRadiationUp",&SEMIC::LongwaveRadiationUp)
		.def("RunEnergyBalance",&SEMIC::RunEnergyBalance)
		.def("RunMassBalance",&SEMIC::RunMassBalance)
		.def("RunEnergyAndMassBalance",&SEMIC::RunEnergyAndMassBalance);
}
