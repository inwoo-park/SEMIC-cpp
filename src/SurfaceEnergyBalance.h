#include<iostream>
#include<vector>
#include<math.h>

const double PI=3.141592653589793238462643;
const double t0=237.15;
const double sigma=5.67; /* Stefan-Bolzmann conant */
const double eps=0.62197; /* ratio of the molar weight of water vapor to the molar weight of dry air */

using namespace std;

class SEMIC{ /* {{{ */
	public:
		/* Initialize variable {{{ */
		int alb_scheme;     /* Albedo scheme: 0: Slater, 1: ISBA, Denby, 2: Alex, 3: None */

		int nx;     /* number of grid points */		
		int n_ksub; /* number of sub-daily time steps */

		double csh;   /* Sensible heat exchange coefficient */
		double cap;   /* air specific heat capacity */

		double ceff; /* Surface specific heat capacity of snow/ice [J K-1 m-2] */
		vector<double> albi; // Background albedo (bare ice) [no unit]
		vector<double> albl; // Background albedo (bare land) [no unit]

		double hcrit; /* critical snow height for which grid cell is 50% snow covered */
		double rcrit; /* Critical snow height for which refreezing fraction is 50% */

		double tmin;  /* Minimum temperature for which albedo decline becomes effective ('slater') [unit: K] */
		double tmax;  /* Maximum temperature for which albedo decline becomes effective ('slater') [unit: K] */

		vector<double> t2m;    /* 2-m air temperature [unit: K] */
		vector<double> tsurf;  /* surface temperature [unit: K] */
		vector<double> hsnow;  /* Snow pack height (water equivalent) [unit: m] */
		vector<double> hice;   /* ice thickness (water equivalent) */
		vector<double> melt;   /* potential surface melt [m/s] */
		vector<double> amp;    /* Temperature Amplitude [unit: K] */
		vector<double> smb;
		vector<double> smb_ice;
		vector<double> smb_snow;
		vector<double> runoff; /* potential surface runoff [m/s] */

		vector<double> rhoa;   /* air density */
		vector<double> wind;   /* Wind speed value [m/s] */

		vector<double> shf;    /* sensibale heat flux [W m-2] */
		vector<double> qmr;    /* heat flux melting/refreezing */
		/* }}} */

		void 	Initialize(int nx);
		void 	Display();
		void 	SensibleHeatFlux();
		void 	LatentHeatFlux(double *sp);
		double 	SaturateWaterVaporP(double temperature);
		void  	RunEnergyBalance();
		void    Run();
}; /* }}} */
