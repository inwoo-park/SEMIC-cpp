#ifndef SURFACEENERGYBALANCE_H
#define SURFACEENERGYBALANCE_H

#include<iostream>
#include<vector>
#include<cmath>
#include<cassert>
#include"SemicParameters.h"

/* Initialize global constant value */
const double PI=3.141592653589793238462643;
const double T0=237.15;
const double SIGMA=5.67e-8; /* Stefan-Bolzmann conant */
const double EPSIL=0.62197; /* ratio of the molar weight of water vapor to the molar weight of dry air */

const double CLS 	= 2.83e+6; /* latent heat of sublimation (J kg-1) */
const double CLM 	= 3.3e+5; 	/* latent heat of melting (J kg-1)*/
const double CLV 	= 2.5e+6; 	/* latent heat of condensation (J kg-1) */
const double CAP 	= 1000.;  	/* specific heat capacity of air (J kg-1 K-1) */
const double RHOW 	= 1000.;  	/* density of water (kg m-3) */
const double HSMAX 	= 5; 		/* maximum snow height (m) */

using namespace std;

/* Type definition */
typedef vector<double> DoubleVector;
typedef vector<int> IntVector;

class SEMIC{ /* {{{ */
	public:
		/* Initialize variable {{{ */
		int alb_scheme;     /* Albedo scheme: 0: Slater, 1: ISBA, Denby, 2: Alex, 3: None */

		int nx;     /* number of grid points */		
		int n_ksub; /* number of sub-daily time steps */

		SemicParameters *Param;
		SemicConstants 	*Const;

		DoubleVector t2m;    /* 2-m air temperature [unit: K] */
		DoubleVector tsurf;  /* surface temperature [unit: K] */
		DoubleVector hsnow;  /* Snow pack height (water equivalent) [unit: m] */
		DoubleVector hice;   /* ice thickness (water equivalent) */
		DoubleVector alb;  /* grid-averaged albedo [no unit] */
		DoubleVector alb_snow; /* snow albedo [not unit] */
		DoubleVector melt;   /* potential surface melt [m/s] */
		DoubleVector melted_snow;   /* actual melted snow [m s-1] */
		DoubleVector melted_ice;    /* actual melted ice [m s-1]*/
		DoubleVector refr;          /* Refreezing [m s-1]*/
		
		DoubleVector smb; 			/* Surface amss balance [m s-1] */
		DoubleVector acc; 			/* surface accumulation [m s-1] */
		DoubleVector lhf; 			/* Latent heat flux [W m-2]*/
		DoubleVector smb_ice;
		DoubleVector smb_snow;
		DoubleVector runoff; /* potential surface runoff [m/s] */

		DoubleVector shf;    /* sensibale heat flux [W m-2] */
		DoubleVector qmr;    /* heat flux melting/refreezing */
		DoubleVector qmr_res;  /* Residual heat flux from melting/refreezing (at end of time step) [W m-2]*/
		
		DoubleVector amp;    /* Temperature Amplitude [unit: K] */
		DoubleVector sf;     /* snowfall [w.e. m s-1] */
		DoubleVector rf;     /* rainfall [w.e. m s-1] */
		DoubleVector sp;     /* Surface pressure [Pa] */
		DoubleVector lwd;    /* Long-wave radiation downward direction [W m-2] */
		DoubleVector swd;    /* Short-wave radiation downward direction [W m-2] */
		DoubleVector wind;   /* Surface wind speed [m s-1] */
		DoubleVector rhoa;   /* Air density [kg m-3] */
		DoubleVector qq;     /* air specific humidity [kg kg-1] */
		IntVector mask;      /* */

		/* auxilary variables */
		DoubleVector lwup; /* upward longwave radiation [W m-2] */
		DoubleVector subl;
		DoubleVector evap;

		/* }}} */

		SEMIC(void);
	    ~SEMIC(void);

		void 	Initialize(int nx);
		void 	Display();
		void 	SensibleHeatFlux(SemicParameters *Param, SemicConstants *Const);
		void 	LatentHeatFlux(SemicParameters Param, SemicConstants Const, DoubleVector sp, DoubleVector wind, DoubleVector &evap, DoubleVector &subl, DoubleVector &lhf);
		double 	SaturateWaterVaporP(double temperature);
		
		void    LongwaveRadiationUp();
		void    TestReturnVector(vector<double> &tmp);

		void DiurnalCycle(DoubleVector &tmean_input, DoubleVector &above, DoubleVector &below);

		/* Albedo schemes */
		double Albedo_Slater(double tsurf, double tmin, double tmax, double alb_smax, double alb_smin);
		double Albedo_Denby(double melt, double alb_smax, double alb_smin, double mcrit);
		double Albedo_ISBA(double alb, double sf, double melt, double tstic, double tau, double tau_a, double tau_f, double w_crit, double mcrit, double alb_smin, double alb_smax);

		double 	ew_sat(double t);
		double 	ei_sat(double t);

		/* Go Solve!*/
		void  	RunEnergyBalance();
		void 	RunMassBalance();
		void 	RunEnergyAndMassBalance();
}; /* }}} */

#endif
