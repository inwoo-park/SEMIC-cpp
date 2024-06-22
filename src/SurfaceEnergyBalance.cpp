#include "omp.h"
#include "SemicParameters.h"
#include "SurfaceEnergyBalance.h"

SEMIC::SEMIC(void){ /*{{{*/
	/* nothing to do. */
	this->Param = new SemicParameters();
	this->Const = new SemicConstants();

	this->verbose = true;

	this->InitializeParameters();
} /*}}}*/
SEMIC::~SEMIC(void){ /*{{{*/
	delete Param;
	delete Const;
} /*}}}*/

void SEMIC::Initialize(int nx){ /* {{{ */
		/*
		* Initialize SEMIC class.
		*/
		this->nx         = nx;
		this->alb_scheme = 0;

		//vector<double> zeros(nx,0);

		this->t2m    = DoubleVector(nx,0.0);
		this->tsurf  = DoubleVector(nx,0.0);
		this->shf    = DoubleVector(nx,0.0);

		this->wind = DoubleVector(nx,0.0); 
		this->rhoa = DoubleVector(nx,0.0);
		this->qq   = DoubleVector(nx,0.0);

		this->lwup = DoubleVector(nx,0.0);

		this->alb      = DoubleVector(nx,0.0);
		this->alb_snow = DoubleVector(nx,0.0);

		this->smb      = DoubleVector(nx,0.0);
		this->smb_snow = DoubleVector(nx,0.0);
		this->smb_ice  = DoubleVector(nx,0.0);

		this->acc		= DoubleVector(nx,0.0);

		this->melt        = DoubleVector(nx,0.0);
		this->melted_snow = DoubleVector(nx,0.0);
		this->melted_ice  = DoubleVector(nx,0.0);

		this->runoff      = DoubleVector(nx,0.0);

		this->refr = DoubleVector(nx,0.0);

		this->hsnow = DoubleVector(nx,0.0);
		this->hice  = DoubleVector(nx,0.0);
		
		this->subl = DoubleVector(nx,0.0);
		this->evap = DoubleVector(nx,0.0);
		this->lhf  = DoubleVector(nx,0.0);

		this->qmr    = DoubleVector(nx,0.0);
		this->qmr_res= DoubleVector(nx,0.0);
} /* }}} */

void SEMIC::InitializeParameters(void){ /*{{{*/
	/* Initialize parameters */
	this->n_ksub = 3;
	this->Param->ceff = 2e+6;
	this->Param->tstic = 86400.;

	this->Param->albi = 0.07;
	this->Param->albl = 0.15;
	this->Param->tmin = -999;
	this->Param->tmax = 273.15;
	this->Param->hcrit = 0.028;
	this->Param->rcrit = 0.79;
	this->Param->amp = DoubleVector(nx,3.5);
	this->alb_scheme = 0;
	this->Param->tau_a = 0.008;
	this->Param->tau_f = 0.24;
	this->Param->w_crit = 0.15;
	this->Param->mcrit = 6e-8;
} /*}}}*/

void SEMIC::Display(){ /* {{{ */
	cout << this->nx << "\n";
	for (int i=0; i<this->smb.size(); i++){
		cout << this->smb[i] << " " << this->t2m[i] << "\n";
	}
} /* }}} */

void SEMIC::SensibleHeatFlux(SemicParameters *Param, SemicConstants *Const){ /* {{{ */
	/* semic-f90: sensible_heat_flux

	*/
	if (this->verbose) cout << "calculate sensible heat flux\n";

	double csh=Param->csh; /* sensible heat exchange coefficient */
	double cap=Const->cap; /* air specific heat capacity */

	for (int i=0; i < this->nx; i++)
		this->shf[i] = csh * cap * this->rhoa[i] * this->wind[i] * (this->tsurf[i] - this->t2m[i]);
} /* }}} */

void SEMIC::LatentHeatFlux(SemicParameters *Param, SemicConstants *Const){/*{{{*/
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
	int i, nx = this->nx;
	DoubleVector shumidity_sat(nx,0); /* Specific humidity */
	DoubleVector sat_vaporP(nx,0);    /* Saturated water vapor pressure*/

	double clh=Param->clh;
	double cls=Const->cls;
	double clv=Const->clv;
	
	/* Check input data size */
	assert(this->wind.size() == nx);
	assert(this->rhoa.size() == nx);
	assert(this->sp.size() == nx);
	assert(this->evap.size() == nx);
	assert(this->subl.size() == nx);
	assert(this->lhf.size() == nx);

	if (this->verbose) cout << "calculate latent heat flux\n";

	for (i=0; i < nx; i++){
		if (this->tsurf[i] < T0){
			sat_vaporP[i] = this->ei_sat(this->tsurf[i]);
			/* specific humidity at surface (assumed to be saturated) is */
			shumidity_sat[i] = sat_vaporP[i] * EPSIL /(sat_vaporP[i]*(EPSIL-1.0)+ this->sp[i]);
			
			this->subl[i] = clh*this->wind[i]*this->rhoa[i]*(shumidity_sat[i] - this->qq[i]);

			this->lhf[i] = this->subl[i]*cls;
		}
		else{
			sat_vaporP[i] = this->ew_sat(this->tsurf[i]);
			/* Evaporation, condenstation */
			/* Specific humidity at surface (assumed to be saturated) */
			shumidity_sat[i] = sat_vaporP[i] * EPSIL /(sat_vaporP[i]*(EPSIL-1.0) + this->sp[i]);
			this->evap[i] = clh*this->wind[i]*this->rhoa[i]*(shumidity_sat[i] - qq[i]);
			this->lhf[i] = this->evap[i]*clv;
		}
	}
} /*}}}*/

void SEMIC::LongwaveRadiationUp(){ /*{{{*/
	/* Calculate upward long-wave radiation with Stefan-Boltzman law

		lwup = sigma * T^4
	*/
	int nx=this->nx;
	
	if (this->verbose) cout << "calculate upward long-wave radiation.\n";

	for (int i=0; i < nx; i++)
		this->lwup[i] = SIGMA*pow(this->tsurf[i], 4.0);
} /*}}}*/

double SEMIC::SaturateWaterVaporP(double temperature) { /* {{{ */
	/* Saturation water vapor pressure over ice 
		
		temperature - temperature (unit: K)

		Output
		fsta - saturated water vapor pressure (unit: P)
	*/
	double fsat; /* saturatoin water vapor */
	fsat = 611.2 * exp(22.46 * (temperature-T0)/(272.62+temperature-T0));

	return fsat;
} /* }}} */

void SEMIC::TestReturnVector(vector<double> &tmp){ /*{{{*/
	for (int i=0; i<tmp.size(); i++){
		tmp[i] = (double)i;
	}
}/*}}}*/

void SEMIC::DiurnalCycle(double tmean, double amp, double &above, double &below){ /*{{{*/
	/* semic-f90: diurnal_cycle
	*/
	double tmp1, tmp2; /* temporal variable */
	
	/* Check consistency*/
	if (abs(tmean/amp) < 1.0){
		tmp1 = acos(tmean/amp);
		tmp2 = sqrt(1 - pow(tmean,2.0)/pow(amp,2.0));
	}

	if (tmean + amp < 0.0){
		below = tmean;
		above = 0.0;
	}else{
		above = tmean;
		below = 0.0;
		if (abs(tmean) < amp){
			/* dt = 2. * x1 */
			above = (-tmean*tmp1 + amp*tmp2 + PI*tmean)/(PI - tmp1);
			/* dt = x2 - x1 */
			below = (tmean*tmp1 - amp*tmp2)/tmp1;
		}
	}
} /*}}}*/

/* Should implent Slater Albedo scheme */
double SEMIC::Albedo_Slater(double tsurf, double tmin, double tmax, double alb_smax, double alb_smin){/*{{{*/
	/* Calculate the albedo with Slater's method

	alb_smax - maximum snow albedo (unit: none)
	alb_smin - minimum snow albedo (unit: none)
	tmin - minimum temperature (unit: K )
	tmax - maximum temperature (unit: K)

	Returns.
	alb - return value
	*/
	double alb;
	double tm = 0.0; /* initialize tm temperature */
	double f; 

	/* flexible factor ensures continuous polynomial */
	f = 1.0/(T0 - tmin);
	if ((tsurf >= tmin) && (tsurf <=tmax)){
		tm = f*(tsurf - tmin);
	}
	if (tsurf > T0){
		tm = 1.0;
	}
	/*
	Krapp et al. (2017)'s comment.
	tm = dmax1(dmin1((tsurf-tmin)/(t0-tmin),1.0_dp),0.0_dp)
	In contrast to the formulation in their paper, I summed up alpha_nir
	and alpha_nir immediately (fewer parameters: alb_smax and alb_smin).
	*/
	alb = alb_smax - (alb_smax - alb_smin)*pow(tm,3.0);
	return alb;
}/*}}}*/
double SEMIC::Albedo_Denby(double melt, double alb_smax, double alb_smin, double mcrit) {/*{{{*/
	/* Calculate Albedo with Denby model 

	Reference
	* 
	*/
	double alb;
	alb = alb_smin + (alb_smax - alb_smin) * exp(-melt/mcrit);
	return alb;
}/*}}}*/
double SEMIC::Albedo_ISBA(double alb, double sf, double melt, double tstic, double tau, double tau_a, double tau_f, double w_crit, double mcrit, double alb_smin, double alb_smax){/*{{{*/
	/* Calcaulate snow albedo with ISBA method
	
	alb - input initial snow albedo
	melt - melting rate
	*/
	double alb_dry, alb_wet, alb_new;
	double w_alb;

	/* where no melting occurs, albedo decreases linearly */
	alb_dry = alb - tau_a*tstic/tau;
	/* where melting occurs, albedo decreases exponentially */
	alb_wet = (alb - alb_smin) * exp(-tau_f*tstic/tau) + alb_smin;
	alb_new = sf*tstic/(w_crit/RHOW) * (alb_smax - alb_smin);

	/* Dry/wet-averaged albedo */
	w_alb = 0.;
	if (melt >0.)
		w_alb = 1 - melt/mcrit;
	w_alb = min(1.0, max(w_alb, 0.));

	alb = (1.0 - w_alb) * alb_dry + w_alb * alb_wet + alb_new;
	alb = min(alb_smax, max(alb, alb_smin));
	return alb;
}/*}}}*/

double SEMIC::ew_sat(double t){/*{{{*/
	/* Saturation water vapor prepssure over water

	t - temperature (unit: K)
	fsat - saturation water vapor (unit: Pa)
	*/
	double fsat;
	fsat = 611.2 * exp(17.62 * (t-T0) /(243.12+t-T0));
	return fsat;
}/*}}}*/
double SEMIC::ei_sat(double t){/*{{{*/
	/* Saturation water vapor prepssure over ice.

	t - temperature (unit: K)
	fsat - saturation water vapor (unit: Pa)
	*/
	double fsat;
	fsat = 611.2 * exp(22.46 * (t-T0) /(272.62+t-T0));
	return fsat;
}/*}}}*/

void SEMIC::RunEnergyBalance() { /* {{{ */
	/* Calculate surface energy balance */
	int i;
	int nx = this->nx; /* size of element*/
	
	/* initialize auxilary variables */
	DoubleVector qsb(nx,0);
	assert(this->mask.size() == nx);
	assert(this->t2m.size() == nx);
	assert(this->shf.size() == nx);
	assert(this->lhf.size() == nx);

	if (this->verbose) cout << "   RunEnergyBalance\n";
	
	/* 1. Calculate the sensible heat flux */
	this->SensibleHeatFlux(this->Param, this->Const);

	/* 2. Calculate the latent heat flux */
	fill(this->subl.begin(), this->subl.end(), 0.0);
	fill(this->evap.begin(), this->evap.end(), 0.0);
	fill(this->lhf.begin(), this->lhf.end(), 0.0);
	this->LatentHeatFlux(this->Param, this->Const);

	#pragma omp parallel for
	for (i=0; i<nx; i++)
		this->subl[i] = this->subl[i]/RHOW;

	/* 3. Surface physics: long-wave radiation */
	assert(this->lwup.size() == nx);
	fill(this->lwup.begin(), this->lwup.end(), 0.0);
	this->LongwaveRadiationUp();

	/* 4. Calculate surface energy balance of incoming and outgoing surface fluxes (W m-2) */
	if (this->verbose) cout << "   step4: calculate surface energy balance\n";
	#pragma omp parallel for
	for (i=0; i<nx; i++){
		qsb[i] = (1. - this->alb[i])*this->swd[i] + this->lwd[i] - this->lwup[i] - this->shf[i] - this->lhf[i] - this->qmr_res[i];
	}

	/* 5. Update surface temperature acoording to surface energy balancec */
	if (this->verbose) cout << "   step5: update surface temperature\n";
	assert(this->qmr.size() == nx);
	assert(this->tsurf.size() == nx);
	#pragma omp parallel for
	for (i = 0; i < nx; i++){
		this->qmr[i] = 0.0; /* set qmr as zeros. */
		this->tsurf[i] = this->tsurf[i] + qsb[i]*this->Param->tsticsub / this->Param->ceff;
		/* store residual energy for subfreezing tsurf over ice and thick snow cover in qmr for later use in mass balance */
		if (((this->mask[i] == 2) || (this->hsnow[i] > 0.0)) && (this->tsurf[i] > T0)){
			this->qmr[i] = (this->tsurf[i]-T0) * this->Param->ceff / this->Param->tsticsub;
			this->tsurf[i] = T0;
		}
	}

	/* 6. Update 2-m air temperature over ice sheet */
	if (this->verbose) cout << "   step6: update 2-air temperature\n";
	#pragma omp parallel for
	for (i = 0; i < nx; i++){
		if ((this->mask[i] == 2) || (this->hsnow[i] > 0)){
			this->t2m[i] = this->t2m[i] + (this->shf[i] + this->lhf[i]) * this->Param->tsticsub / this->Param->ceff;
		}
	}

	if (this->verbose) cout << "   Finalize clear memory\n";	
	qsb.clear();
} /* }}} */

void SEMIC::RunMassBalance(){/*{{{*/
	/* Calculate mass balance
	*/
	int i;
	int nx = this->nx;
	DoubleVector above(nx,0.0), below(nx,0.0);
	DoubleVector qmelt(nx,0.0), qcold(nx,0.0);
	DoubleVector snow_to_ice_input(nx);
	DoubleVector refrozen_rain(nx);
	DoubleVector refrozen_snow(nx);

	assert(this->melt.size() == nx);
	assert(this->melted_snow.size() == nx);
	assert(this->melted_ice.size() == nx);
	assert(this->hsnow.size() == nx);
	assert(this->acc.size() == nx);
	assert(this->smb_snow.size() == nx);
	assert(this->sf.size() == nx);
	assert(this->subl.size() == nx);
	assert(this->smb.size() == nx);
	assert(this->smb_snow.size() == nx);
	assert(this->smb_ice.size() == nx);

	if (this->verbose) cout << "   RunMassBalance\n";

	#pragma omp parallel for
	for (i=0; i<nx; i++){
		/* 1. Calculate above/below freezing temperature for a given mean temeprature*/
		if (this->verbose) cout << "   step1: calculate above/below freezing temperature\n";
		this->DiurnalCycle(this->tsurf[i] - T0 + this->qmr[i]/this->Param->ceff*this->Param->tstic,
					this->Param->amp[i], above[i], below[i]);

		this->qmr[i] = 0.;
	
		if (this->mask[i] >= 1){
			/* 2. Calculate melt energy where temperature exceeds freezing */	
			qmelt[i] = max(0.0, above[i] * this->Param->ceff / this->Param->tstic);
			/* 3. Calcaulate "cold" content */
			qcold[i] = max(0.0, abs(below[i]) * this->Param->ceff / this->Param->tstic);
		}
		else{
			qmelt[i] = 0.;
			qcold[i] = 0.;
		}

		/* 4. Ablation: melt (water m s-1); potential melt resulting from available melt energy*/
		if (this->verbose) cout << "   step4: calculate ablation.\n";

		/* potential melt */
		this->melt[i] = qmelt[i]/(RHOW * CLM);
		/* seperate potential melt into actual melt of snow and ice */

		this->melted_snow[i] = min(this->melt[i], this->hsnow[i]/this->Param->tstic);
		
		this->melted_ice[i] = this->melt[i] - this->melted_snow[i];

		if (this->mask[i] == 2){
			/* actual melt is sum of melted snow and ice (melted snow over land) */
			this->melt[i] = this->melted_snow[i] + this->melted_ice[i];
		}
		else{
			this->melt[i] = this->melted_snow[i];
			this->melted_ice[i] = 0;
		}

		/* 5. Refreezing (m s-1) as fraction of melt (increase with snow height) */
		if (this->verbose) cout << "   step5: refreezing.\n";
		double f_rz = this->Param->rcrit; /* get freezing parameter*/

		this->refr[i] = qcold[i] / (RHOW * CLM);
		/* potential refreezing */
		refrozen_rain[i] = min(this->refr[i], this->rf[i]);
		/* potential refreezing snow */
		refrozen_snow[i] = min(this->refr[i] - refrozen_rain[i], 0.);
		/* actual refreezing snow*/
		refrozen_snow[i] = min(refrozen_snow[i], this->melted_snow[i]);
		/* actual refreezing */
		refrozen_rain[i] = f_rz * min(this->hsnow[i]/this->Param->tstic, refrozen_rain[i]);
		refrozen_snow[i] = f_rz * min(this->hsnow[i]/this->Param->tstic, refrozen_snow[i]);
		this->refr[i] = refrozen_rain[i] + refrozen_snow[i];
		/* energy released during refreezing that has not been used
		   is subtracted from residual energy
		*/
		this->qmr[i] = this->qmr[i] - (1.0 - f_rz) * this->refr[i] * RHOW *CLM;

		/* 6. Runoff */
		if (this->verbose) cout << "   step6: runoff\n";
		this->runoff[i] = this->melt[i] + this->rf[i] - refrozen_rain[i];

		/* 7. Accumulation: sum of all incoming solid water (just diagnostic here) */
		if (this->verbose) cout << "   step7: accumulation\n";
		this->acc[i] = this->sf[i] - this->subl[i] + this->refr[i];

		/* 8. Surface mass balancee of snow */
		if (this->verbose) cout << "   step8: surface mass balance of snow\n";
		this->smb_snow[i] = this->sf[i] - this->subl[i] - this->melted_snow[i] + refrozen_snow[i];

		if (this->mask[i] == 0){
			this->hsnow[i] = 0.;
		}
		else{
			/* 9. Update snow height*/
			this->hsnow[i] = max(0.0, this->hsnow[i] + this->smb_snow[i] * this->Param->tstic);
		}

		/* 10. Relax snow height to maximum (e.g., 5m) */
		if (this->verbose) cout << "   step10: relax snow height\n";
		double snow_to_ice;

		snow_to_ice = max(0.0, this->hsnow[i] - HSMAX);
		/* save temporal data set.*/
		snow_to_ice_input[i] = snow_to_ice;

		this->hsnow[i] = this->hsnow[i] - snow_to_ice;
		this->smb_ice[i] = (snow_to_ice)/(this->Param->tstic) - this->melted_ice[i] + refrozen_rain[i];

		/* use to force ice sheet model */
		/* update new ice budge remove or add ice */
		this->hice[i] = this->hice[i] + this->smb_ice[i] * this->Param->tstic;

		/* 11. Total surface mass balance */
		if (this->verbose) cout << "   step11: total surface mass balance\n";
		snow_to_ice = snow_to_ice_input[i];
		if (this->mask[i] == 2){
			this->smb[i] = this->smb_snow[i] + this->smb_ice[i] - snow_to_ice/this->Param->tstic;
		}
		else{
			this->smb[i] = this->smb_snow[i] + max(0.0, this->smb_ice[i] - snow_to_ice/this->Param->tstic);
		}

		/* 12. Update snow albedo*/
		if (this->verbose) cout << "   step12: update snow albedo\n";
		double f_alb;
		double albi=this->Param->albi;
		double albl=this->Param->albl;
		f_alb = 1 - exp(-this->hsnow[i]/(this->Param->hcrit + EPSILON));
		if (this->alb_scheme == 0){
			this->alb_snow[i] = this->Param->alb_smax;
		}
		else if (this->alb_scheme == 1){
			/* Slater's albedo scheme */
			this->alb_snow[i] = this->Albedo_Slater(this->tsurf[i], this->Param->tmin,
			this->Param->tmax, this->Param->alb_smax, this->Param->alb_smin);
		}
		else if (this->alb_scheme == 2){
			/* Denby albedo scheme */
			this->alb_snow[i] = this->Albedo_Denby(this->melt[i], this->Param->alb_smax, this->Param->alb_smin, this->Param->mcrit);
		}
		else if (this->alb_scheme == 3){
			/* ISBA albedo scheme */
			this->alb_snow[i] = this->Albedo_ISBA(this->alb_snow[i], this->sf[i], this->melt[i], this->Param->tstic, this->Param->tstic, this->Param->tau_f, this->Param->tau_f, this->Param->w_crit, this->Param->mcrit, this->Param->alb_smax, this->Param->alb_smax);
		}
		else{
			cerr << "ERROR: Given alb_scheme is not available\n";
		}

		switch (this->mask[i])
		{
		case 2: /* for ice */
			this->alb[i] = albi + f_alb*(this->alb_snow[i] - albi);
			break;
		case 1: /* for land */
			this->alb[i] = albl + f_alb*(this->alb_snow[i] - albl);
			break;
		case 0: /* for ocean */
			this->alb[i] = 0.06;
			break;
		default:
			cerr << "ERROR: Given mask is not supported.";
			break;
		}

		/* Store residual energy */
		if (this->mask[i]==0) /* for ocean mask*/
			this->qmr_res[i] = 0.0;
		else
			this->qmr_res[i] = this->qmr[i];
	}

	/* Clear memory */
	below.clear();
	above.clear();
	qmelt.clear();
	qcold.clear();

	refrozen_rain.clear();
}/*}}}*/

void SEMIC::RunEnergyAndMassBalance(){ /*{{{*/
	/* Run 
	*/
	int i, nx=this->nx;

	/* sub time stepping */
	this->Param->tsticsub = this->Param->tstic/this->n_ksub;

	//if (this->verbose){
	//	cout << "n_ksub value = " << this->n_ksub << endl;;
	//	cout << "tsticsub value = " << this->Param->tsticsub << endl;
	//}

	if (this->verbose) cout << "Run Energy Balance\n";
	for (i=0; i < this->n_ksub; i++){
		this->RunEnergyBalance();
	}

	if (this->verbose) cout << "Run Mass Balance\n";
	this->RunMassBalance();
} /* }}} */
