#ifndef SEMICPARAMETERS_H
#define SEMICPARAMETERS_H

#include "iostream"
#include "vector"

/* Type definition */
typedef std::vector<double> DoubleVector;
typedef std::vector<int> IntVector;

using namespace std;

class SemicParameters{ /* {{{ */
	public:
		double ceff; /* Surface specific heat capacity of snow/ice [J K-1 m-2] */
		double albi; /* Background albedo (bare ice) [no unit] */
		double albl; /* Background albedo (bare land) [no unit] */
		double alb_smax; /* Maximum snow albedo (fresh snow) [no unit] */
		double alb_smin; /* Minimum snow albedo (fresh snow) [no unit] */
		
		double hcrit; /* Critical snow height for which grid cell is 50% snow covered */
		double rcrit; /* Critical snow height for which refreezing fraction is 50% */
		DoubleVector amp;   /* Amplitude of diuranl cycle [unit: K] */
		double csh;   /* Sensible heat exchange coefficient */
		double clh;   /* Latent heat exchange coefficient [no unit] */
		
		double tmin;  /* Minimum temperature for which albedo decline becomes effective ('slater') [unit: K] */
		double tmax;  /* Maximum temperature for which albedo decline becomes effective ('slater') [unit: K] */

		double tstic; 		/* Time step [unit: s]*/
		double tsticsub; 	/* Sub-time step [unit: s] */

		double tau_a; /* Dry-albedo decline for "isba" albedo scheme [1 day-1] */
		double tau_f; /* Wet albedo declinte for "isba" albedo scheme [1 day=1] */

		double w_crit; /* critical liquid water content for "isba" albedo scheme. [unit: kg m-2]*/

		double mcrit;  /* critical melt rate for "isba" and "denby" albedo scheme [unit: m s-1]*/
		double tmid;   /* parameter for "alex" albedo parameterization [unit: K]*/
		// double clv;   /* latent heat of sublimation */
		// double cls;   /* latent heat of vaporization */
		// double cap;   /* air specific heat capacity */
}; /* }}} */

class SemicConstants{ /* {{{ */
public:
	double cls=2.83e+6; /* latent heat of sublimation [unit: J kg-1] */
	double clm=3.3e+5;  /* latent heat of melting [unit: J kg-1] */
	double clv=2.5e+6;  /* latent heat of condensation [unit: J kg-1] */
	double cap=1e+3;  	/* specific heat capacity of air [unit: J kg-1 K-1]*/
	double rhow=1e+3; 	/* density of water [kg m-3]*/
	double hsmax=5; 	/* maximum snow height [m] */
	// double epsil
	
	void Display(void){
		cout << "cls    " << this->cls  << " latent heat of sublimation (unit: J kg-)" << endl;
		cout << "clm    " << this->clm  << " latent heat of melting (unit: J kg-)" << endl;
		cout << "clv    " << this->clv  << " latent heat of condensation (unit: J kg-)" << endl;
		cout << "cap    " << this->cap  << " specific heat capacity of air (unit: J kg-1 K-1)" << endl;
		cout << "rhow   " << this->rhow << " density of water (unit: kg m-3)" << endl;
	};
}; /* }}}*/
#endif
