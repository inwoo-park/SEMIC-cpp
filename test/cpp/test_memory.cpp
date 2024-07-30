#include "SemicParameters.h"

int main(void){

	int nx = 100;
	int ntime = 1000;

	DoubleMatrix *m = new DoubleMatrix(nx, ntime);
	delete m;

	SemicForcings* f = new SemicForcings(nx, ntime);
	delete f;
	return 0;
}
