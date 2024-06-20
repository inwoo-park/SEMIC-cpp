// #include<iostream>
#include "SurfaceEnergyBalance.h"

void testFunction(double &k){
	k = -1000.0;
}

int main(void){

	/* Initialize SEMIC module */
	SEMIC *semic = new SEMIC();
	int nx=10;
	semic->Initialize(nx);

	vector<double> smb(nx,1);
	vector<double> t2m(nx,273.15);
	vector<double> tsurf(nx,273.15);
	vector<double> lwup(nx);

	for (int i=0; i < nx; i++){
		tsurf[i] = tsurf[i] - (double)i;
	}

	semic->smb = smb;
	semic->t2m = t2m;
	semic->tsurf = tsurf;

	// semic.Display();
	// semic.LongwaveRadiationUp(lwup);
	for (int i=0; i<nx; i++){
		cout << semic->lwup[i] << " " << semic->tsurf[i] << endl;
	}

	semic->LongwaveRadiationUp();
	for (int i=0; i<nx; i++){
		cout << semic->lwup[i] << endl;
	}

	cout << "Show Param values.\n"; 
	cout << semic->Param->alb_smax << endl;
	
	/* set parameters */
	semic->Param->alb_smax = 0.8;
	cout << semic->Param->alb_smax << endl;
	
	delete semic;

	testFunction(tsurf[0]);
	for (int i=0;i<nx; i++)
		cout << tsurf[i] << endl;

	return 0;
}

