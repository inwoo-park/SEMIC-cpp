#include "SemicParameters.h"
#include "SurfaceEnergyBalance.h"

int main(void){

	int nx = 20000;
	int ntime = 7000;

	DoubleMatrix *m = new DoubleMatrix(nx, ntime);
	delete m;

	SemicForcings* f = new SemicForcings(nx, ntime);
	delete f;

	SemicResult* r = new SemicResult(nx, ntime);
	delete r;

	SEMIC* semic = new SEMIC();
	semic->Initialize(nx);

	delete semic;
	return 0;
}
