#include "SurfaceEnergyBalance.h"
#include "SemicParameters.h"

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

		this->wind = zeros; 
		this->rhoa = zeros;
		this->qq   = zeros;
} /* }}} */

void SEMIC::Display(){ /* {{{ */
	cout << this->nx << "\n";
	for (int i=0; i<this->smb.size(); i++){
		cout << this->smb[i] << " " << this->t2m[i] << "\n";
	}
} /* }}} */

void SEMIC::SensibleHeatFlux(SemicParameters Param, SemicConstants Const){ /* {{{ */
	/* semic-f90: sensible_heat_flux

	*/
	double csh=Param.csh; /* sensible heat exchange coefficient */
	double cap=Const.cap; /* air specific heat capacity */
	// csh=Param.csh;
	// cap=Param.cap;

	DoubleVector shf(this->nx, 0);

	for (int i=0; i < this->nx; i++)
		shf[i] = csh * cap * this->rhoa[i] * this->wind[i] * (this->tsurf[i] - this->t2m[i]);

	// Update shf results
	this->shf = shf;
} /* }}} */

void SEMIC::LatentHeatFlux(SemicParameters Param, SemicConstants Const, DoubleVector sp, DoubleVector wind, DoubleVector &evap, DoubleVector &subl, DoubleVector &lhf){ /* {{{ */
	/* semic-f90: latent_heat_flux
	Inputs
	----------
	sp   - surface pressure (unit: Pa).
	wind - wind speed (unit: m s-1)

	Ouputs
	------
	evap    - evapotranspiration
	subl    - sublimation
	lhf     - latent heat flux
	*/
	vector<double> shumidity_sat(this->nx,0); /* Specific humidity */
	vector<double> sat_vaporP(this->nx,0); /* Saturated water vapor pressure*/

	double clh=Param.clh;
	double clv=Const.clv;
	double cls=Const.cls;
	
	/* Check input data size */
	assert(sp.size() == this->nx);
	assert(evap.size() == this->nx);
	assert(subl.size() == this->nx);
	assert(lhf.size() == this->nx);

	for (int i=0; i<this->nx; i++){
		if (this->tsurf[i] < t0){
			sat_vaporP[i] = this->SaturateWaterVaporP(this->tsurf[i]);

			/* Calculate specific humidity*/
			shumidity_sat[i] = sat_vaporP[i]*eps/(sat_vaporP[i]*(eps-1.0)+ sp[i]);
			
			subl[i] = clh*this->wind[i]*(shumidity_sat[i] - this->qq[i]);

			lhf[i] = subl[i]*cls;
		}
		else{
			sat_vaporP[i] = this->SaturateWaterVaporP(this->tsurf[i]);
			/* Evaporation, condenstation */
			/* Specific humidity at surface (assumed to be saturated) */
			shumidity_sat[i] = sat_vaporP[i] * eps /(sat_vaporP[i]*(eps) + sp[i]);
			evap[i] = clh*this->wind[i]*this->rhoa[i]*(shumidity_sat[i] - qq[i]);
			lhf[i] = evap[i]*clv;
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

void SEMIC::LongwaveRadiationUp(vector<double> &lwup){
	/* Calculate upward long-wave radiation with Stefan-Boltzman law
	Ouputs
	------
	lwup - upward longwave radiation.
	*/
	int nx=this->nx;
	
	for (int i=0; i < nx; i++)
		lwup[i] = sigma*pow(this->tsurf[i], 4);
}

void SEMIC::TestReturnVector(vector<double> &tmp){
	for (int i=0; i<tmp.size(); i++){
		tmp[i] = (double)i;
	}
}

void SEMIC::DiurnalCycle(SemicParameters Param, DoubleVector &above, DoubleVector &below){
	/* semic-f90: diurnal_cycle
	*/
	double amp = Param.amp; /* temperature amplitude*/
	double tmean; /* surface temperature */
	double tmp1, tmp2; /* temporal variable */
	
	/* Check consistency*/
	assert(above.size() == this->nx);
	assert(below.size() == this->nx);
	
	for (int i=0; i<this->nx; i++){
		tmean = this->tsurf[i]; /* get surface temperature */

		if (abs(tmean/amp) < 1.0){
			tmp1 = acos(tmean/amp);
			tmp2 = sqrt(1 - tmean);
		}

		if (tmean + amp < 0.0){
			below[i] = tmean;
			above[i] = 0.0;
		}else{
			above[i] = tmean;
			below[i] = 0.0;
			if (abs(tmean) < amp){
				/* dt = 2. * x1 */
				above[i] = (-tmean*tmp1 + amp*tmp2 + PI*tmean)/(PI - tmp1);
				/* dt = x2 - x1 */
				below[i] = (tmean*tmp1 - amp*tmp2)/tmp1;
			}
		}
		
	}
}

void SEMIC::RunEnergyBalance() { /* {{{ */
	int nx = this->nx; /* size of element*/
	
	// SemicParameters Param=this->Param;

	DoubleVector lwup(nx,0);

	/* 1. Calculate the sensible heat flux */
	this->SensibleHeatFlux(this->Param, this->Const);

	/* 2. Calculate the latent heat flux */
	// this->LatentHeatFlux();

	/* 3. Surface physics: long-wave radiation */
	this->LongwaveRadiationUp(lwup);

	/* 4. Calculate surface energy balance of incoming and outgoing surfaec flux (W m-2)*/
	// for (i=0; i<nx; i++){

	// }
} /* }}} */

void SEMIC::Run(){ /* {{{ */
	/* Run main function.
		*/

	/* Run energy balance */


	/* Run masss balance */
	} /* }}} */
