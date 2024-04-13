#include<iostream>
#include"SurfaceMassBalance.h"

int main(void){

	/* Initialize SEMIC module */
	SEMIC semic;
	int nx=10;
	semic.Initialize(nx);

	vector<double> smb(nx,1);
	vector<double> t2m(nx,273.15);

	semic.smb = smb;
	semic.t2m = t2m;

	semic.Display();

	return 0;
}

