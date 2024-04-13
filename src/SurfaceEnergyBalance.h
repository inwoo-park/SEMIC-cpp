#ifndef SURFACEENERGYBALANCE_H
#define SURFACEENERGYBALANCE_H

#include<iostream>
#include<vector>
#include<math.h>
#include<cassert>
#include"SemicParameters.h"

const double PI=3.141592653589793238462643;
const double t0=237.15;
const double sigma=5.67; /* Stefan-Bolzmann conant */
const double eps=0.62197; /* ratio of the molar weight of water vapor to the molar weight of dry air */

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

		SemicParameters 	Param;
		SemicConstants 		Const;

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
		/* }}} */

		void 	Initialize(int nx);
		void 	Display();
		void SensibleHeatFlux(SemicParameters Param, SemicConstants Const);
		void LatentHeatFlux(SemicParameters Param, SemicConstants Const, DoubleVector sp, DoubleVector wind, DoubleVector &evap, DoubleVector &subl, DoubleVector &lhf);
		double 	SaturateWaterVaporP(double temperature);
		
		void    LongwaveRadiationUp(vector<double> &lwup);
		void    TestReturnVector(vector<double> &tmp);

		void DiurnalCycle(SemicParameters Param, DoubleVector &above, DoubleVector &below);

		/* Go Solve!*/
		void  	RunEnergyBalance();
		void    Run();

		
}; /* }}} */

#endif
