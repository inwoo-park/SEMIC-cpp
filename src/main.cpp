// #include<iostream>
#include "SurfaceEnergyBalance.h"

int main(void){

	/* Initialize SEMIC module */
	SEMIC semic;
	int nx=10;
	semic.Initialize(nx);

	vector<double> smb(nx,1);
	vector<double> t2m(nx,273.15);
	vector<double> tsurf(nx,273.15);
	vector<double> lwup(nx);

	semic.smb = smb;
	semic.t2m = t2m;
	semic.tsurf = tsurf;

	// semic.Display();
	// semic.LongwaveRadiationUp(lwup);
	semic.TestReturnVector(lwup);
	for (int i=0; i<nx; i++){
		cout << lwup[i] << endl;
	}

	return 0;
}

