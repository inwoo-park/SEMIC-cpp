#include "SemicParameters.h"
#include "SurfaceEnergyBalance.h"

void test_class_pointer(){ /* {{{ */
	int nx = 10000;
	int ntime = 1000;

	DoubleMatrix *m = new DoubleMatrix(nx, ntime);
	delete m;

	SemicForcings* f = new SemicForcings(nx, ntime);
	delete f;

	SemicResult* r = new SemicResult(nx, ntime);
	delete r;

	SEMIC* semic = new SEMIC();
	semic->Initialize(nx);

	delete semic;
} /* }}} */

void test_DoubleMatrix(){ /* {{{ */
	DoubleMatrix *m = new DoubleMatrix();

	cout << "Initialize m matrix\n";
	cout << "nx, ny = " << m->nrow << "," << m->ncol << endl;

	vector<vector<double>> a(20, vector<double>(30,0.0));
	vector<vector<double>> b;

	/* reshape array */
	m->set_value(a);
	b = m->get_value();

	cout << "Set value, m matrix\n";
	cout << "nx, ny = " << m->nrow << "," << m->ncol << endl;

	cout << "Information for b matrix\n";
	cout << "nx, ny = " << b.size() << "," << b[0].size() << endl;

	delete m;
} /* }}} */

int main(void){
	//test_class_pointer();
	test_DoubleMatrix();
	return 0;
}
