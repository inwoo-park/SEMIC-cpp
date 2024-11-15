#include "omp.h"
#include "SemicParameters.h"
#include "SurfaceEnergyBalance.h"

/* Check Elapsed time! */
#include <chrono>

SEMIC::SEMIC(void){ /*{{{*/
	/* nothing to do. */
	this->Param = new SemicParameters();
	this->Const = new SemicConstants();
   this->Result = new SemicResult();

	this->verbose = true;

	this->num_threads = 1; /* use single cpus */
	this->SetOpenmpThreads(1); /* default setting with single cpu usage.*/
} /*}}}*/
SEMIC::~SEMIC(void){ /*{{{*/
	delete this->Param;
	delete this->Const;
   delete this->Result;
} /*}}}*/

void SEMIC::Initialize(int nx){ /* {{{ */
		/*
		* Initialize SEMIC class.
		*/
		this->nx         = nx;

		double initValue = 0.0;
		this->InitializeParameters();
		this->alb_scheme = 0;

		this->t2m    = DoubleVector(nx,initValue);
		this->tsurf  = DoubleVector(nx,initValue);
		this->shf    = DoubleVector(nx,initValue);

		this->wind = DoubleVector(nx,initValue); 
		this->rhoa = DoubleVector(nx,initValue);
		this->qq   = DoubleVector(nx,initValue);

		this->lwup = DoubleVector(nx,initValue);

		this->alb      = DoubleVector(nx,initValue);
		this->alb_snow = DoubleVector(nx,initValue);

		this->smb      = DoubleVector(nx,initValue);
		this->smb_snow = DoubleVector(nx,initValue);
		this->smb_ice  = DoubleVector(nx,initValue);

		this->acc		= DoubleVector(nx,initValue);

		this->melt        = DoubleVector(nx,initValue);
		this->melted_snow = DoubleVector(nx,initValue);
		this->melted_ice  = DoubleVector(nx,initValue);

		this->runoff      = DoubleVector(nx,initValue);

		this->refr = DoubleVector(nx,initValue);

		this->hsnow = DoubleVector(nx,initValue);
		this->hice  = DoubleVector(nx,initValue);
		
		this->subl = DoubleVector(nx,initValue);
		this->evap = DoubleVector(nx,initValue);
		this->lhf  = DoubleVector(nx,initValue);

		this->qmr    = DoubleVector(nx,initValue);
		this->qmr_res= DoubleVector(nx,initValue);
} /* }}} */

void SEMIC::InitializeParameters(void){ /*{{{*/
	/* first check nx value! */
	if (this->nx == 0){
		throw runtime_error("ERROR: SEMIC->nx is not defined!");
	}
	double init_amp = 3.5;

	/* Initialize parameters */
	this->n_ksub = 3;
	this->Param->tstic = 86400.;

	this->Param->albi = 0.07;
	this->Param->albl = 0.15;

    /* use default albedo integration with Oerlemans et al. (1998) method

    * Oerlemans, J., & Knap, W. H. (1998). A 1 year record of global radiation and albedo in the ablation zone of Morteratschgletscher, Switzerland. Journal of Glaciology, 44(147), 231–238. https://doi.org/10.3189/S0022143000002574
    */
    this->Param->alb_scheme_sum = 0;

	/* for slater's albedo scheme 
	Use tmin with -10 C noted in Slater's 1998 paper.
	*/
	this->Param->tmin = 273.15-10;
	this->Param->tmax = 273.15;

	this->Param->hcrit = 0.028;
	this->Param->rcrit = 0.79;
	this->Param->amp = DoubleVector(this->nx, init_amp);
	this->alb_scheme = 1; /* Use slater's albedo scheme */
	this->Param->tau_a = 0.008;
	this->Param->tau_f = 0.24;
	this->Param->w_crit = 0.15;
	this->Param->mcrit = 6e-8;
} /*}}}*/

void SEMIC::Display(){ /* {{{ */
	cout << this->nx << "\n";
	for (int i=0; i<(int)this->smb.size(); i++){
		cout << this->smb[i] << " " << this->t2m[i] << "\n";
	}
} /* }}} */

void SEMIC::SensibleHeatFlux(SemicParameters *Param, SemicConstants *Const){ /* {{{ */
	/* semic-f90: sensible_heat_flux
	Calculate sensible heat flux.

	*/
	// if (this->verbose) cout << "calculate sensible heat flux\n";

	double csh=Param->csh; /* sensible heat exchange coefficient */
	double cap=Const->cap; /* air specific heat capacity */

	// #pragma omp parallel for schedule(dynamic)
	for (int i=0; i < this->nx; i++)
		this->shf[i] = csh * cap * this->rhoa[i] * this->wind[i] * (this->tsurf[i] - this->t2m[i]);
} /* }}} */

double SEMIC::SensibleHeatFlux(SemicParameters *Param, SemicConstants *Const, double rhoa, double wind, double tsurf, double t2m) {/* {{{ */
	/* semic-f90: sensible_heat_flux
	Calcaulate sensible heat flux with given variable

	shf - return value: sensible heat flux
	*/
	double csh=Param->csh; /* sensible heat exchange coefficient */
	double cap=Const->cap; /* air specific heat capacity */
	double shf=0.0; /* sensible heat flux */
	
	/* sensible heat flux */
	shf = csh * cap * rhoa * wind * (tsurf - t2m);
	return shf;
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

	// #pragma omp parallel for schedule(dynamic)
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

void SEMIC::LatentHeatFlux(SemicParameters *Param, SemicConstants *Const, double rhoa, double wind, double tsurf, double sp, double qq, double &evap, double &subl, double &lhf){/*{{{*/
	/* semic-f90: latent_heat_flux
	Inputs
	----------
	sp   - surface pressure (unit: Pa).
	wind - wind speed (unit: m s-1)

	Ouputs (returns)
	------
	evap    - evapotranspiration
	subl    - sublimation
	lhf     - latent heat flux
	*/
	double shumidity_sat=0.0; /* Specific humidity */
	double sat_vaporP=0.0; /* Saturated water vapor pressure */

	/* From Parameters */
	double clh=Param->clh;

	/* From Constants */
	double cls=Const->cls; 
	double clv=Const->clv;
	
	/* Initialize variable */
	subl = 0.0;
	lhf  = 0.0;
	evap = 0.0;

	if (tsurf < T0){
		sat_vaporP = this->ei_sat(tsurf);
		/* specific humidity at surface (assumed to be saturated) is */
		shumidity_sat = sat_vaporP * EPSIL /(sat_vaporP*(EPSIL-1.0)+ sp);
		
		subl= clh*wind*rhoa*(shumidity_sat - qq);
		lhf = subl*cls;
	}
	else{
		sat_vaporP = this->ew_sat(tsurf);
		/* Evaporation, condenstation */
		/* Specific humidity at surface (assumed to be saturated) */
		shumidity_sat = sat_vaporP * EPSIL /(sat_vaporP*(EPSIL-1.0) + sp);
		evap = clh*wind*rhoa*(shumidity_sat - qq);
		lhf  = evap*clv;
	}
} /*}}}*/

void SEMIC::LongwaveRadiationUp(){ /*{{{*/
	/* Calculate upward long-wave radiation with Stefan-Boltzman law

		lwup = sigma * T^4
	*/
	int nx=this->nx;
	
	if (this->verbose) cout << "calculate upward long-wave radiation.\n";

	for (int i=0; i < nx; i++)
		this->LongwaveRadiationUp(tsurf[i], this->lwup[i]);
} /*}}}*/

void SEMIC::LongwaveRadiationUp(double tsurf, double &lwup){ /*{{{*/
	/* Calculate upward long-wave radiation with Stefan-Boltzman law

		lwup = sigma * T^4

	Inputs
	------
	tsurf - surface temperature (unit: K)

	Outputs (returns)
	-----------------
	lwup - upward long wave radiation (unit: W m-2)
	*/
	// if (this->verbose) cout << "calculate upward long-wave radiation.\n";

	lwup = SIGMA*pow(tsurf, 4.0);
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
	for (int i=0; i<(int)tmp.size(); i++){
		tmp[i] = (double)i;
	}
}/*}}}*/

void SEMIC::DiurnalCycle(double tmean, double amp, double &above, double &below){ /*{{{*/
	/* semic-f90: diurnal_cycle

	Out
	---
	above: double - above temperataure (unit: K)
	below: double - below temperature (unit: K)
	*/
	double tmp1=0.0;
	double tmp2=0.0; /* temporal variable */
	
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
	tmin - minimum temperature (unit: K)
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
	/*
	@brief Calcaulate snow albedo with ISBA method. It requires previous step of snow albedo.
	
	@param alb - input initial snow albedo (unit: m s-1)
	@param melt - melting rate (unit :m s-1)
	@param sf  - snowfall flux (unit: m s-1)
	@param melt - surface melting value (unit: m s-1)
	@param tau_a, tau_f - parameters with value of 0.008 and 0.24, respectively.
	@w_crit - 

	@return alb - surface albedo (unit: -)
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
double SEMIC::Albedo_Zeitz(double melt, double alb_smax, double alb_smin, double mcrit) { /* {{{ */
	/* Calculate snow albedo with Zeitz et al. (2021), Garbe et al. (2023)
	 
		References
		* Zeitz, M., Reese, R., Beckmann, J., Krebs-Kanzow, U., and Winkelmann, R.: Impact of the melt–albedo feedback on the future evolution of the Greenland Ice Sheet with PISM-dEBM-simple, The Cryosphere, 15, 5739–5764, https://doi.org/10.5194/tc-15-5739-2021, 2021.
		* Garbe, J., Zeitz, M., Krebs-Kanzow, U., and Winkelmann, R.: The evolution of future Antarctic surface melt using PISM-dEBM-simple, The Cryosphere, 17, 4571–4599, https://doi.org/10.5194/tc-17-4571-2023, 2023.
	 */
	
	double alb;
	alb = alb_smax - mcrit * melt;
	alb = max(alb, alb_smin);

	return alb;
}/* }}} */

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

void SEMIC::InitializeSemicResult(int nx, int ntime){ /* {{{ */
    /* Initialize new Result class {{{ */
	// if (!this->Result)
	//	 this->Result = new SemicResult(); /* Initialize emtpy array */

	 cout << "InitializeSemicResult: Use output array with set_dimension.\n";

    if (this->Result->iscontain(this->output_request, "smb"))
		  cout << "InitializeSemicResult: set smb array\n";
        //this->Result->smb = new DoubleMatrix(nx,ntime);
        this->Result->smb->set_dimension(nx,ntime);
    if (this->Result->iscontain(this->output_request, "melt"))
        //this->Result->melt = new DoubleMatrix(nx,ntime);
        this->Result->melt->set_dimension(nx,ntime);
    if (this->Result->iscontain(this->output_request, "alb"))
        //this->Result->alb= new DoubleMatrix(nx,ntime);
        this->Result->alb->set_dimension(nx,ntime);
    if (this->Result->iscontain(this->output_request, "tsurf"))
        //this->Result->tsurf= new DoubleMatrix(nx,ntime);
        this->Result->tsurf->set_dimension(nx,ntime);
    if (this->Result->iscontain(this->output_request, "hsnow"))
        //this->Result->hsnow= new DoubleMatrix(nx,ntime);
        this->Result->hsnow->set_dimension(nx,ntime);
	 /* }}} */
} /* }}} */

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

	/* Initialize default value */
	fill(this->subl.begin(), this->subl.end(), 0.0);
	fill(this->evap.begin(), this->evap.end(), 0.0);
	fill(this->lhf.begin(), this->lhf.end(), 0.0);

	if (this->verbose)
		cout << "   RunEnergyBalance\n";

	assert(this->lwup.size() == nx);
	assert(this->qmr.size() == nx);
	assert(this->tsurf.size() == nx);

	fill(this->lwup.begin(), this->lwup.end(), 0.0);
	
	int thread_id = 0;
	
	#pragma omp parallel private(i, thread_id) shared(nx, qsb) /* Initialize Openmp Parallel*/
	{
	 thread_id = omp_get_thread_num();

	//#pragma omp master
	if ((this->verbose) & (thread_id == 0))
		cout << "   step1: Calculate the sensible heat flux\n";
	
	#pragma omp for //schedule(dynamic)
	for (int i=0; i<nx; i++){

		this->shf[i] = this->SensibleHeatFlux(this->Param, this->Const, this->rhoa[i], this->wind[i], this->tsurf[i], this->t2m[i]);
	}

	if (this->verbose && thread_id == 0)
		cout << "   step2: Calculate the latent heat flux\n";

	#pragma omp for //schedule(dynamic)	
	for (int i=0; i<nx; i++){
		this->LatentHeatFlux(this->Param, this->Const, this->rhoa[i], this->wind[i], this->tsurf[i], this->sp[i], this->qq[i], this->evap[i], this->subl[i], this->lhf[i]);

		this->subl[i] = this->subl[i]/RHOW;
	}
	
	//#pragma omp master
	if (this->verbose && thread_id == 0)
		cout << "   step3: Surface physics: long-wave radiation\n";
	
	#pragma omp for
	for (int i=0; i<nx; i++)
		this->LongwaveRadiationUp(this->tsurf[i], this->lwup[i]);
	
	/* 4. Calculate surface energy balance of incoming and outgoing surface fluxes (W m-2) */
	//#pragma omp master
	if ((this->verbose) & (thread_id == 0))
		cout << "   step4: calculate surface energy balance\n";

	#pragma omp for
	for (int i=0; i<nx; i++){
		qsb[i] = (1. - this->alb[i])*this->swd[i] + this->lwd[i] - this->lwup[i] - this->shf[i] - this->lhf[i] - this->qmr_res[i];
	}

	/* 5. Update surface temperature acoording to surface energy balancec */
	//#pragma omp master
	if ((this->verbose) & (thread_id == 9))
		cout << "   step5: update surface temperature\n";

	#pragma omp for 
	for (int i=0; i<nx; i++){

		this->qmr[i] = 0.0; /* set qmr as zeros. */
		this->tsurf[i] = this->tsurf[i] + qsb[i]*this->Param->tsticsub / this->Param->ceff;
		/* store residual energy for subfreezing tsurf over ice and thick snow cover in qmr for later use in mass balance */
		if (((this->mask[i] == 2) || (this->hsnow[i] > 0.0)) && (this->tsurf[i] > T0)){
			this->qmr[i] = (this->tsurf[i]-T0) * this->Param->ceff / this->Param->tsticsub;
			this->tsurf[i] = T0;
		}
	}

	/* 6. Update 2-m air temperature over ice sheet */
	#pragma omp master
	{
	if (this->verbose)
		cout << "   step6: update 2-air temperature\n";
	}

	#pragma omp for //schedule(dynamic)
	for (int i=0; i<nx; i++){
		if ((this->mask[i] == 2) || (this->hsnow[i] > 0.0)){
			this->t2m[i] = this->t2m[i] + (this->shf[i] + this->lhf[i]) * this->Param->tsticsub / this->Param->ceff;
		}
	}

	} /* End of Openmp Parallel */

	if (this->verbose) cout << "   Finalize clear memory\n";	
	qsb.clear();
} /* }}} */

void SEMIC::RunMassBalance(){/*{{{*/
	/* Calculate mass balance
	*/
	int i;
	int nx = this->nx;
	double f_rz; 
	DoubleVector above(nx,0.0), below(nx,0.0);
	DoubleVector qmelt(nx,0.0), qcold(nx,0.0);
	DoubleVector snow_to_ice_input(nx, 0.0);
	DoubleVector refrozen_rain(nx, 0.0);
	DoubleVector refrozen_snow(nx, 0.0);

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

	int thread_id = 0;

	if (this->verbose) cout << "   RunMassBalance\n";

    /* start omp parallel*/
	#pragma omp parallel private(i, f_rz, thread_id) shared(nx, snow_to_ice_input, qmelt, qcold, above, below, refrozen_rain, refrozen_snow)
	{
	 thread_id = omp_get_thread_num(); /* get master thread */

	#pragma omp for //schedule(dynamic)
	for (i=0; i<nx; i++){
		/* 1. Calculate above/below freezing temperature for a given mean temeprature*/
		if (i==0 && this->verbose){
			cout << "   step1: calculate above/below freezing temperature\n";
			cout << "      tsurf= " << this->tsurf[i] << endl;
			cout << "      qmr  = " << this->qmr[i] << endl;
			cout << "      Tamp = " << this->Param->amp[i] << endl;
		}
		this->DiurnalCycle(this->tsurf[i] - T0 + this->qmr[i]/this->Param->ceff * this->Param->tstic,
					this->Param->amp[i], above[i], below[i]);
		
		if (i==0 && this->verbose){
			cout << "      above temprature " << above[i] << endl;
			cout << "      below temprature " << below[i] << endl;
		}

		if (i==0 && this->verbose)
			cout << "   -- Set qmr = 0.\n";
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
		if (i==0 && this->verbose) cout << "   step4: calculate ablation.\n";

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
		if (i==0 && this->verbose) cout << "   step5: refreezing.\n";
		f_rz = this->Param->rcrit; /* get freezing parameter*/

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
		if (i==0 && this->verbose) cout << "   step6: runoff\n";
		this->runoff[i] = this->melt[i] + this->rf[i] - refrozen_rain[i];

		/* 7. Accumulation: sum of all incoming solid water (just diagnostic here) */
		if (i==0 && this->verbose) cout << "   step7: accumulation\n";
		this->acc[i] = this->sf[i] - this->subl[i] + this->refr[i];

		/* 8. Surface mass balancee of snow */
		if (i==0 && this->verbose)
			cout << "   step8: surface mass balance of snow\n";
		this->smb_snow[i] = this->sf[i] - this->subl[i] - this->melted_snow[i] + refrozen_snow[i];

		if (this->mask[i] == 0){
			this->hsnow[i] = 0.;
		}
		else{
			/* 9. Update snow height*/
			this->hsnow[i] = max(0.0, this->hsnow[i] + this->smb_snow[i] * this->Param->tstic);
		}

		/* 10. Relax snow height to maximum (e.g., 5m) */
		if (i==0 && this->verbose) cout << "   step10: relax snow height\n";
		double snow_to_ice;

		snow_to_ice = max(0.0, this->hsnow[i] - HSMAX);
		/* save temporal data set.*/
		//snow_to_ice_input[i] = snow_to_ice;

		this->hsnow[i] = this->hsnow[i] - snow_to_ice;
		this->smb_ice[i] = (snow_to_ice)/(this->Param->tstic) - this->melted_ice[i] + refrozen_rain[i];

		/* use to force ice sheet model */
		/* update new ice budge remove or add ice */
		this->hice[i] = this->hice[i] + this->smb_ice[i] * this->Param->tstic;

		/* 11. Total surface mass balance */
		if (i==0 && this->verbose)
			cout << "   step11: total surface mass balance\n";
		//snow_to_ice = snow_to_ice_input[i];
		if (this->mask[i] == 2){
			this->smb[i] = this->smb_snow[i] + this->smb_ice[i] - snow_to_ice/this->Param->tstic;
		}
		else{
			this->smb[i] = this->smb_snow[i] + max(0.0, this->smb_ice[i] - snow_to_ice/this->Param->tstic);
		}

		/* 12. Update snow albedo*/
		if (i==0 && this->verbose)
			cout << "   step12: update snow albedo\n";
		
		double f_alb=0.0;
		double albi=this->Param->albi;
		double albl=this->Param->albl;

		/* Calculate */
		if (this->Param->alb_scheme_sum == 0)
			f_alb = 1 - exp(-this->hsnow[i]/(this->Param->hcrit + EPSILON));
		else if (this->Param->alb_scheme_sum == 1){
			/* Eq. (11) in Napoly et al. (2020)
			
			Napoly, A., Boone, A., & Welfringer, T. (2020). ISBA-MEB (SURFEX v8.1): Model snow evaluation for local-scale forest sites. Geoscientific Model Development, 13(12), 6523–6545. https://doi.org/10.5194/gmd-13-6523-2020
			*/
			f_alb = 1 - this->hsnow[i]/(this->Param->hcrit + EPSILON);
			/* preventing negative value! */
			f_alb = max(0.0, min(1.0, f_alb)); /* 0<= f_alb <= 1 */
		}
		else{
			cerr << "ERROR: Given Param->alb_scheme_sum (=" << this->Param->alb_scheme_sum << ") value is not available\n";
		}

		/* Now, calcaulate snow albedo*/
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
		else if (this->alb_scheme == 4){
			/* Zetiz albedo scheme */
			this->alb_snow[i] = this->Albedo_Zeitz(this->melt[i], this->Param->alb_smax, this->Param->alb_smin, this->Param->mcrit);
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

	} /* end of pragma omp*/

	/* Clear memory */
	below.clear();
	above.clear();
	qmelt.clear();
	qcold.clear();

	refrozen_rain.clear();
}/*}}}*/

void SEMIC::RunEnergyAndMassBalance(){ /*{{{*/
	/* Run both energy and mass balance
	*/
	int i;

	/* sub time stepping */
	this->Param->tsticsub = this->Param->tstic/this->n_ksub;

	/* Check # of threads */
	/*
	#pragma omp parallel
	{
		#pragma omp master
		cout << "Number of threads = " << omp_get_num_threads() << endl;
	}
	*/

	if (this->verbose) cout << "Run Energy Balance\n";
	for (i=0; i < this->n_ksub; i++){
		this->RunEnergyBalance();
	}

	if (this->verbose) cout << "Run Mass Balance\n";
	this->RunMassBalance();
} /* }}} */

void SEMIC::RunEnergyAndMassBalance(SemicForcings* Forcings, int nloop){ /* {{{ */
    /* Run SEMIC with SemicForcing class. */

	int i, j;
	int nx=this->nx;
	int ntime=Forcings->ntime;

    /* Variables for checking elapsed time */
    chrono::system_clock::time_point tstart;
    chrono::duration<double> dt1, dt2;

    /* Check consistency*/
    assert(this->nx == Forcings->nx);

    /* Initialize new Result class {{{ */
    if (this->Result->iscontain(this->output_request, "smb")){
#ifdef HAVE_DEBUG
            cout << "SemicResults: Set SMB!\n";
#endif
        this->Result->smb->set_dimension(nx, ntime);
    }
    if (this->Result->iscontain(this->output_request, "melt")){
#ifdef HAVE_DEBUG
            cout << "SemicResults: Set Melt!\n";
#endif
        this->Result->melt->set_dimension(nx,ntime);
    }
    if (this->Result->iscontain(this->output_request, "alb")){
#ifdef HAVE_DEBUG
            cout << "SemicResults: Set Alb!\n";
#endif
        this->Result->alb->set_dimension(nx,ntime);
    }
    if (this->Result->iscontain(this->output_request, "smb_snow"))
        this->Result->smb_snow->set_dimension(nx,ntime);
    if (this->Result->iscontain(this->output_request, "smb_ice"))
        this->Result->smb_ice->set_dimension(nx,ntime);
    if (this->Result->iscontain(this->output_request, "alb_snow"))
        this->Result->alb_snow->set_dimension(nx,ntime);
    if (this->Result->iscontain(this->output_request, "tsurf"))
        this->Result->tsurf->set_dimension(nx,ntime);
    if (this->Result->iscontain(this->output_request, "hsnow"))
        this->Result->hsnow->set_dimension(nx,ntime);
    if (this->Result->iscontain(this->output_request, "hice"))
        this->Result->hice->set_dimension(nx,ntime);
    if (this->Result->iscontain(this->output_request, "subl"))
        this->Result->subl->set_dimension(nx,ntime);
    if (this->Result->iscontain(this->output_request, "evap"))
        this->Result->evap->set_dimension(nx,ntime);

	 /* }}} */

	 /* sub-time stepping */
	 this->Param->tsticsub = this->Param->tstic/this->n_ksub;

	//  DoubleVector sf, rf, sp;
	//  DoubleVector lwd, swd, wind, rhoa, t2m, qq;

    for (int _nloop=0; _nloop<nloop; _nloop++){
        for (j = 0; j<ntime; j++){
            if (this->verbose)
                cout << "Time step " << j+1 << "/" << ntime << endl;
            
            /* Update forcing varibales */
            if (this->verbose)
                cout << "   Assign variable" << endl;
            this->rf   = Forcings->rf->get_column_vector(j);
            this->sf   = Forcings->sf->get_column_vector(j);
            this->sp   = Forcings->sp->get_column_vector(j);
            
            this->lwd  = Forcings->lwd->get_column_vector(j);
            this->swd  = Forcings->swd->get_column_vector(j);
            this->wind = Forcings->wind->get_column_vector(j);
            this->rhoa = Forcings->rhoa->get_column_vector(j);
            this->t2m  = Forcings->t2m->get_column_vector(j);
            this->qq   = Forcings->qq->get_column_vector(j);

            // #ifdef HAVE_DEBUG
            //     if (j == 0 && _nloop == 0)
            //     cout << "nrow and ncol = (" << Forcings->t2m->nrow << " " << Forcings->t2m->ncol << ")\n";
            // #endif
            
            // #ifdef HAVE_DEBUG
            //     /* check forcing varible */
            //     if (_nloop == 0){
            //         cout << "ntime = " << j << ", t2m[" << 0 << "] = " << this->t2m[0] <<  " Forcings->t2m[0] = " << Forcings->t2m->value[0][j] << endl;
            //     }
            // #endif

            /* Now, use forcings variables! */
            if (this->verbose) cout << "Run Energy Balance\n";
            tstart = chrono::system_clock::now();
            for (i=0; i < this->n_ksub; i++){
                this->RunEnergyBalance();
            }
            dt1 = chrono::system_clock::now() - tstart;

            if (this->verbose) cout << "Run Mass Balance\n";
            tstart = chrono::system_clock::now();
            this->RunMassBalance();
            dt2 = chrono::system_clock::now() - tstart;

            if (this->verbose){
                cout << "Elapsed time\n";
                cout << "-- Energy balance = " << dt1.count() << endl;
                cout << "-- Mass balance   = " << dt2.count() << endl;
            }

            /* Return output value */
            if (_nloop == nloop-1){
#ifdef HAVE_DEBUG
                if (j == 0)
                    cout << "Export Results!\n";
#endif
                for (i=0; i<nx; i++){
                    if (this->Result->iscontain(this->output_request, "smb"))
                        this->Result->smb->value[i][j]   = this->smb[i];
                    if (this->Result->iscontain(this->output_request, "smb_snow"))
                        this->Result->smb_snow->value[i][j]   = this->smb_snow[i];
                    if (this->Result->iscontain(this->output_request, "smb_ice"))
                        this->Result->smb_ice->value[i][j]   = this->smb_ice[i];
                    if (this->Result->iscontain(this->output_request, "melt"))
                        this->Result->melt->value[i][j]  = this->melt[i];
                    if (this->Result->iscontain(this->output_request, "tsurf"))
                        this->Result->tsurf->value[i][j] = this->tsurf[i];
                    if (this->Result->iscontain(this->output_request, "alb"))
                        this->Result->alb->value[i][j]   = this->alb[i];
                    if (this->Result->iscontain(this->output_request, "alb_snow"))
                        this->Result->alb_snow->value[i][j]   = this->alb_snow[i];
                    if (this->Result->iscontain(this->output_request, "hsnow"))
                        this->Result->hsnow->value[i][j]   = this->hsnow[i];
                    if (this->Result->iscontain(this->output_request, "hice"))
                        this->Result->hice->value[i][j]   = this->hice[i];
                    if (this->Result->iscontain(this->output_request, "subl"))
                        this->Result->subl->value[i][j]   = this->subl[i];
                    if (this->Result->iscontain(this->output_request, "evap"))
                        this->Result->evap->value[i][j]   = this->evap[i];
                }
            }
        }
    }
} /* }}} */

void SEMIC::RunEnergyAndMassBalance(SemicForcings* Forcings){ /* {{{ */
    /* Only run one-time */
    this->RunEnergyAndMassBalance(Forcings, 1);
} /* }}} */

void SEMIC::SetOpenmpThreads(void){ /* {{{ */ 
	omp_set_num_threads(this->num_threads);
} /* }}} */

void SEMIC::SetOpenmpThreads(int ncpus){ /* {{{ */
	omp_set_num_threads(ncpus);
} /* }}} */

int SEMIC::GetOpenmpThreads(void){ /* {{{ */ 
	int nthreads = 1;
	#pragma omp parallel
	{
		#pragma omp single
		{
			nthreads = omp_get_num_threads();
		}
	}
	return nthreads;
} /* }}}*/

// void SEMIC::SetOpenmpRuntime(int scheduleType, int chunkSize){ /* {{{ */
// 	/*Initialize schedule type for OpenMP

// 	scheduleType
// 	1: static
// 	2: dynamic
// 	3: guidded
// 	4: auto

// 	chunkSize - chunk size for each thread.
// 	*/
// 	omp_set_schedule(scheduleType, chunkSize);
// } /* }}} */

int SEMIC::GetOpenmpVersion(void){
	/*return OpenMP version.
	*/
	return _OPENMP;
}
