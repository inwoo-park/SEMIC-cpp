#ifndef DATAARRAY_H
#define DATAARRAY_H

#include <vector>
#include <iostream>
#include <stdexcept>

#define USE_TEMPLATE 0
#define _DEBUG_ 0

using namespace std;

class DoubleMatrix{
	public:
		int nrow;
		int ncol;
		double** value;

		/* Initialize constructor*/
		DoubleMatrix(): nrow(0), ncol(0), value(nullptr) {}
		DoubleMatrix(int nrow, int ncol){
			int i,j;
			this->nrow = nrow;
			this->ncol = ncol;

            if (nrow == 0 && ncol == 0){
                /* Initialize null pointer in this->value */
                this->value = nullptr;
            }else{
                AllocateMemory(this->nrow, this->ncol);
			
                for (i=0; i<this->nrow; i++){
                    for (j=0; j<this->ncol; j++){
                        this->value[i][j] = 0.0;
                    }
                }
            }			
		}
		DoubleMatrix(const vector<vector<double>>& values){
			if (!values.size()){
				throw runtime_error("We cannot find dimension in values!");
			}
			
			this->nrow = (int)values.size();
			this->ncol = (int)values[0].size();

			this->AllocateMemory(this->nrow, this->ncol);

			for (int i=0; i<this->nrow; i++){
				for (int j=0; j<this->ncol; j++){
					this->value[i][j] = values[i][j];
				}
			}
		}
		/* Copy constructor */
		DoubleMatrix(const DoubleMatrix& other) : nrow(other.nrow), ncol(other.ncol) {
        	AllocateMemory(this->nrow, this->ncol);
			for (int i = 0; i < nrow; ++i) {
				for (int j = 0; j < ncol; ++j) {
					value[i][j] = other.value[i][j];
				}
			}
    	}

		~DoubleMatrix(){
			DeallocateMemory();
#if _DEBUG_
			std::cout << "Destructor: deallocated memory" << std::endl;
#endif
		}

		void Destroy(){
			DeallocateMemory();
		}

		/* Overload operator() for element access */
		double& operator()(int row, int col) {
			if (row >= nrow || col >= ncol) {
				throw out_of_range("Index out of range");
			}
			return value[row][col];
		}
		const double& operator()(int row, int col) const {
			if (row >= nrow || col >= ncol) {
				throw out_of_range("Index out of range");
			}
			return value[row][col];
		}
		double* operator[](int row){
			if (row < 0 || row >= nrow) {
				throw std::out_of_range("Row index out of range");
			}
			return value[row];
		}
		DoubleMatrix operator+(const DoubleMatrix &B){
			if (this->nrow != B.nrow || this->ncol != B.ncol){
				throw runtime_error("Dimension mismatch");
			}

			DoubleMatrix result = DoubleMatrix(this->nrow, this->ncol);
			for (int i=0; i<this->nrow; i++){
				for (int j=0; j<this->ncol; j++){
					result.value[i][j] = this->value[i][j] + B.value[i][j];
				}
			}
			return result;
		}

		/* Copy assignment operator */
		DoubleMatrix& operator=(const DoubleMatrix& other) {
			if (this == &other) {
				return *this;
			}
			DeallocateMemory();
			nrow = other.nrow;
			ncol = other.ncol;
			AllocateMemory(nrow, ncol);
			for (int i = 0; i < nrow; ++i) {
				for (int j = 0; j < ncol; ++j) {
					value[i][j] = other.value[i][j];
				}
			}
			return *this;
			}

		/* Set values. */
		void set_value(const vector<vector<double>>& values){
			/* Set array's value and matrix size */
			if (this->nrow == 0 && this->ncol == 0){
				this->nrow = (int)values.size();
				this->ncol = (int)values[0].size();

				if (!this->value){
					AllocateMemory(this->nrow, this->ncol);
				}
				for (int i=0; i<this->nrow; i++){
					for (int j=0; j<this->ncol; j++){
						this->value[i][j] = values[i][j];
					}
				}
			} else{
				if (values.size() !=this->nrow || values[0].size() != this->ncol){
					throw runtime_error("Dimension mismatch!");
				}
				if (!this->value) {
					throw std::runtime_error("Memory for matrix not allocated");
				}
				for (int i=0; i<this->nrow; i++){
					for (int j=0; j<this->ncol; j++){
						this->value[i][j] = values[i][j];
					}
				}
			}
		}
		void set_value(const double A){
			int i,j;
			if (!this->value) {
				throw std::runtime_error("Memory for matrix not allocated");
			}

			for (i=0; i<this->nrow; i++){
				for (j=0; j<this->ncol; j++)
					this->value[i][j] = A;
			}
		}

		vector<vector<double>> get_value() const{
			if (!this->value) {
				throw std::runtime_error("Memory for matrix not allocated");
			}

			vector<vector<double>> result(this->nrow, vector<double>(this->ncol));
			for (int i = 0; i<this->nrow; i++){
				for (int j=0; j<this->ncol; j++){
					result[i][j] = value[i][j];
				}
			}
			return result;
		}

		/* To preventing overloading problem, we try to use "template". */
#if USE_TEMPLATE
		template <typename Template>
		Template get_column(int col) const{
			/*Return column vector.
			*/
			if (col > this->ncol || col < 0){
				throw runtime_error("Given col is out of range.");
			}

			if (is_same<Template, double*>::value) {
				double* result = new double[this->nrow];
				for (int i=0; i<this->nrow; ++i){
					result[i] = this->value[i][col];
				}
				return result;
			} else if (is_same<Template, vector<double>>::value) {
				vector<double> result(this->nrow,0.0);
				for (int i=0; i<this->nrow; ++i){
					result[i] = this->value[i][col];
				}
				return result;
			} else {
            	static_assert(always_false<T>::value, "Unsupported return type");
        	}
		}
#else
		double* get_column(int col) const{
			/*Get column vector in "double*" type.
			*/
			if (col > this->ncol || col < 0){
				throw runtime_error("Given col is out of range.");
			}
			double* result = new double[this->nrow];
			for (int i=0; i<this->nrow; ++i){
				result[i] = this->value[i][col];
			}
			return result;	
		}
        vector<double> get_column_vector(int col) const{
			/*Get column vector in "double*" type.
			*/
			if (col > this->ncol || col < 0){
				throw runtime_error("Given col is out of range.");
			}
			vector<double> result(this->nrow);
            for (int i=0; i<this->nrow; ++i){
				result[i] = this->value[i][col];
			}
			return result;	
		}
#endif

		double* get_row(int row) const{
			if (row > this->nrow || row< 0){
				throw runtime_error("Given row is out of range.");
			}

			double* result = new double[this->ncol];
			for (int i=0; i<this->ncol; i++){
				result[i] = this->value[row][i];
			}

			return result;			
		}
		vector<double> get_row_vector(int row) const{
			if (row > this->nrow || row< 0){
				throw runtime_error("Given row is out of range.");
			}

			vector<double> result(this->ncol);
			for (int i=0; i<this->ncol; i++){
				result[i] = this->value[row][i];
			}

			return result;
		}
		
		double sum(){ /* {{{ */
			double gSum = 0.0;
			int i,j;
			for (i=0; i<this->nrow; i++){
				for (j=0; j<this->ncol; j++){
					gSum+=this->value[i][j];
				}
			}
			return gSum;
		} /* }}} */

		double mean(){ /* {{{ */
			return this->sum()/(this->nrow*this->ncol);
		} /* }}} */

		void Display(){ /* {{{ */
			int i, j;
			cout << "Matrix\n";
			for (i=0; i<this->nrow; i++){
				for (j=0; j<this->ncol; j++){
					cout << this->value[i][j] << " ";
				}
				cout << "\n";
			}
		} /* }}} */
private:
	void AllocateMemory(int nrow, int ncol){
		/* Allocate memory*/
		this->value = new double*[nrow];
		for (int i=0; i<this->nrow; i++){
			this->value[i] = new double[this->ncol];
		}
	}

	void DeallocateMemory(){
		int i=0;
		for (i=0; i<this->nrow; i++){
			delete[] value[i];
		}
		delete[] value;
	}

	template <class>
    struct always_false : std::false_type {};
};

#endif