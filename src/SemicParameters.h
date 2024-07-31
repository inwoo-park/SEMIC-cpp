#ifndef SEMICPARAMETERS_H
#define SEMICPARAMETERS_H

#include <iostream>
#include <array> 
#include <vector>
#include <algorithm> /* for using "find" */
// #include "SemicArray.h"
#include "DataArray.h"

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

        SemicParameters(){
            this->amp = DoubleVector(0, 0.0);
        }
        ~SemicParameters(){
            this->amp.clear();
        }
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

class SemicForcings{ /* {{{ */
	public:
		int nx=0; /* numer of grid*/
		int ntime=0; /* number of time step */

		/* set variables */
		DoubleMatrix *sf   = new DoubleMatrix();  /* snowfall [w.e. m s-1] */
		DoubleMatrix *rf   = new DoubleMatrix();  /* rainfall [w.e. m s-1] */
		DoubleMatrix *t2m  = new DoubleMatrix();  /* surface 2-m air temperature [unit: K]*/
		DoubleMatrix *sp   = new DoubleMatrix();  /* Surface pressure [Pa] */
		DoubleMatrix *lwd  = new DoubleMatrix();  /* Long-wave radiation downward direction [W m-2] */
		DoubleMatrix *swd  = new DoubleMatrix();  /* Short-wave radiation downward direction [W m-2] */
		DoubleMatrix *wind = new DoubleMatrix();  /* Surface wind speed [m s-1] */
		DoubleMatrix *rhoa = new DoubleMatrix();  /* Air density [kg m-3] */
		DoubleMatrix *qq   = new DoubleMatrix();  /* air specific humidity [kg kg-1] */	

		/* Initialize constructor and destructor. */
		SemicForcings() {}
		/*
		SemicForcings(): nx(0), ntime(0),
            sf(new DoubleMatrix()),
            rf(new DoubleMatrix()),
            t2m(new DoubleMatrix()),
            sp(new DoubleMatrix()),
            lwd(new DoubleMatrix()),
            swd(new DoubleMatrix()),
            wind(new DoubleMatrix()),
            rhoa(new DoubleMatrix()),
            qq(new DoubleMatrix())            
            {}
		*/

		SemicForcings(int nx, int ntime):
				nx(nx), ntime(ntime),
            sf(new DoubleMatrix(nx, ntime)),
            rf(new DoubleMatrix(nx, ntime)),
            t2m(new DoubleMatrix(nx, ntime)),
            sp(new DoubleMatrix(nx, ntime)),

            lwd(new DoubleMatrix(nx, ntime)),
            swd(new DoubleMatrix(nx, ntime)),
            wind(new DoubleMatrix(nx, ntime)),
            rhoa(new DoubleMatrix(nx, ntime)),
            qq(new DoubleMatrix(nx, ntime)) {}
			
        /* Deconstructor */
		~SemicForcings(){
            delete this->sf;
            delete this->rf;
            delete this->t2m;
            delete this->sp;

            delete this->lwd;
            delete this->swd;
            delete this->wind;
            delete this->rhoa;
            delete this->qq;
        }

		void Initialize(int nx, int ntime){ /* {{{ */
			this->nx = nx;
			this->ntime = ntime;

         this->sf  = new DoubleMatrix(nx, ntime);
         this->rf  = new DoubleMatrix(nx, ntime);
         this->t2m = new DoubleMatrix(nx, ntime);
         this->sp  = new DoubleMatrix(nx, ntime);

         this->lwd = new DoubleMatrix(nx, ntime);
         this->swd = new DoubleMatrix(nx, ntime);

         this->wind = new DoubleMatrix(nx, ntime);
         this->rhoa = new DoubleMatrix(nx, ntime);
         this->qq   = new DoubleMatrix(nx, ntime);
		}; /* }}} */
}; /* }}} */

class SemicResult{ /* {{{ */
    public:
        vector<string> output_request;
        vector<string> output_list = {"smb", "smb_snow","smb_ice",
				"tsurf","melt","alb","hsnow","hice"};
        DoubleMatrix *smb;
        DoubleMatrix *smb_ice;
        DoubleMatrix *smb_snow;
        DoubleMatrix *melt;
        DoubleMatrix *alb;
        DoubleMatrix *alb_snow;
        DoubleMatrix *hsnow;
        DoubleMatrix *hice;

        DoubleMatrix *subl;
        DoubleMatrix *evpl;
        DoubleMatrix *tsurf;

        SemicResult(){ /* {{{ */
            this->smb      = new DoubleMatrix();
            this->smb_ice  = new DoubleMatrix();
            this->smb_snow = new DoubleMatrix();

            this->alb      = new DoubleMatrix();
            this->alb_snow = new DoubleMatrix();

            this->melt  = new DoubleMatrix();
            this->tsurf = new DoubleMatrix();

            this->hsnow = new DoubleMatrix();
            this->hice  = new DoubleMatrix();

            this->subl = new DoubleMatrix();
            this->evpl = new DoubleMatrix();

            /* Initialize output requests */
            this->output_request.push_back("smb");
            this->output_request.push_back("smb_snow");
            this->output_request.push_back("smb_ice");
            this->output_request.push_back("melt");
            this->output_request.push_back("alb");
            this->output_request.push_back("tsurf");
            this->output_request.push_back("hsnow");
            this->output_request.push_back("hice");

            //this->output_list = {"smb", "smb_snow",
            //"smb_ice",
            //"tsurf","melt","alb","hsnow","hice"};
        } /* }}} */

        SemicResult(int nrow, int ncol){ /* {{{ */
            this->smb      = new DoubleMatrix(nrow, ncol);
            this->smb_ice  = new DoubleMatrix(nrow, ncol);
            this->smb_snow = new DoubleMatrix(nrow, ncol);

            this->alb      = new DoubleMatrix(nrow, ncol);
            this->alb_snow = new DoubleMatrix(nrow, ncol);
            
            this->melt     = new DoubleMatrix(nrow, ncol);
            this->tsurf    = new DoubleMatrix(nrow, ncol);
            this->hsnow    = new DoubleMatrix(nrow, ncol);
            this->hice    = new DoubleMatrix(nrow, ncol);

            this->subl = new DoubleMatrix(nrow, ncol);
            this->evpl = new DoubleMatrix(nrow, ncol);
        } /* }}} */

        ~SemicResult(){ /* {{{ */
#ifdef HAVE_DEBUG
            cout << "SemicResult: Destroy memory!\n";
#endif

#ifdef HAVE_DEBUG
			   cout << "Destroy memory: smb, smb_ice, and smb_snow\n";
#endif
            delete this->smb;
            delete this->smb_ice;
            delete this->smb_snow;

#ifdef HAVE_DEBUG
            cout << "Destroy memory: alb\n";
#endif

            delete this->alb;
            delete this->alb_snow;

#ifdef HAVE_DEBUG
			cout << "Destroy memory: melt\n";
#endif
            delete this->melt;
            delete this->tsurf;

            delete this->hsnow;
            delete this->hice;

#ifdef HAVE_DEBUG
			cout << "Destroy memory: subl and evpl\n";
#endif
            delete this->subl;
            delete this->evpl;
        } /* }}} */

        void free_memory(){
#ifdef HAVE_DEBUG
			   cout << "SemicResult: Destroy memory!\n";
#endif

#ifdef HAVE_DEBUG
			   cout << "Destroy memory: smb_ice\n";
#endif
				 delete this->smb_ice;
#ifdef HAVE_DEBUG
			   cout << "Destroy memory: smb_snow\n";
#endif
				 delete this->smb_snow;

#ifdef HAVE_DEBUG
			   cout << "Destroy memory: alb\n";
#endif
				 delete this->alb;

				 delete this->alb_snow;

#ifdef HAVE_DEBUG
			   cout << "Destroy memory: smb_melt\n";
#endif
				 delete this->melt;
				 delete this->tsurf;
		  }

		   /* search if specific string in string array */
		   bool iscontain(vector<string> list, string value){
				return find(list.begin(), list.end(), value) != list.end();
			}

			bool ContainOutput(string value){
				return find(this->output_list.begin(), this->output_list.end(), value) != output_list.end();
			}
}; /* }}} */

#endif
