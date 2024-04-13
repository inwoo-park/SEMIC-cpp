#include "SurfaceEnergyBalance.h"

void SEMIC::Initialize(int nx){ /* {{{ */
		/*
			* Initialize SEMIC class.
			*/
		this->nx         = nx;
		this->alb_scheme = 0;

		vector<double> zeros(nx,0);

		this->t2m   = zeros;
		this->tsurf = zeros;
		this->qmr   = zeros;
		this->shf   = zeros;
		this->smb   = zeros;
} /* }}} */

void SEMIC::Display(){ /* {{{ */
	cout << this->nx << "\n";
	for (int i=0; i<this->smb.size(); i++){
		cout << this->smb[i] << " " << this->t2m[i] << "\n";
	}
} /* }}} */

void SEMIC::SensibleHeatFlux(){ /* {{{ */
	double csh; /* sensible heat exchange coefficient */
	vector<double> shf(this->nx, 0);

	for (int i = 0; i < this->nx; i++)
		shf[i] = csh * cap * this->rhoa[i] * this->wind[i] * (this->tsurf[i] - this->t2m[i]);

	// Update shf results
	this->shf = shf;
} /* }}} */

void SEMIC::LatentHeatFlux(double *sp){ /* {{{ */
	/*
	
	Parameters
	----------
	sp  - surface pressure (unit: Pa).
		*/
	vector<double> shumidity(this->nx); /* Specific humidity */
	vector<double> sat_vaporP(this->nx); /* saturated water vapor pressure*/

	for (int i; i<this->nx; i++){
		if (this->tsurf[i] < t0){
			sat_vaporP[i] = this->SaturateWaterVaporP(this->tsurf[i]);

			/* Calculate specific humidity*/
			shumidity[i] = sat_vaporP[i]*eps/(sat_vaporP[i]*(eps-1.0)+ sp[i]);
		}
		else{

		}
	}
} /* }}} */

double SEMIC::SaturateWaterVaporP(double temperature) { /* {{{ */
	/* Saturation water vapor pressure over ice 
		*
		* temperature - temperature (unit: K)
	*/
	double fsat; /* saturatoin water vapor */
	fsat = 611.2 * exp(22.46 * (temperature-t0)/(272.62+temperature-t0));

	return fsat;
} /* }}} */

void SEMIC::RunEnergyBalance() { /* {{{ */

	/* 1. Calculate the sensible heat flux */
	this->SensibleHeatFlux();

	/* 2. Calculate the latent heat flux */
	this->LatentHeatFlux();
} /* }}} */

void SEMIC::Run(){ /* {{{ */
	/* Run main function.
		*/

	/* Run energy balance */


	/* Run masss balance */
	} /* }}} */
