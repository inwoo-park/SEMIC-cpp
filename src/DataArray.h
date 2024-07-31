#ifndef DATAARRAY_H
#define DATAARRAY_H

#include <vector>
#include <iostream>
#include <stdexcept>

#define USE_TEMPLATE 0
#ifdef HAVE_DEBUG
	#define _DEBUG_ 1
#endif

#ifdef HAVE_PYBIND11
	#include "pybind11/pybind11.h"
	#include "pybind11/stl.h"
	#include "pybind11/numpy.h"
	namespace py = pybind11;
#endif

using namespace std;

class DoubleMatrix{
	public:
		int nrow=0;
		int ncol=0;
		//double** value=nullptr;
		vector<vector<double>> value;

		/* Initialize constructor*/
		DoubleMatrix(void){}  //nrow(0), ncol(0), value(nullptr) {}
		DoubleMatrix(int nrow, int ncol){
			int i,j;
			this->nrow = nrow;
			this->ncol = ncol;

            if (nrow == 0 && ncol == 0){
                /* Initialize null pointer in this->value */
                //this->value = nullptr;
					 /* nothing to do */
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
#ifdef HAVE_PYBIND11
		DoubleMatrix(py::array_t<double, py::array::c_style>& arr) { /* {{{ */
			py::buffer_info info = arr.request();
			auto buf = arr.mutable_unchecked();

			if (info.ndim != 2) {
				throw std::runtime_error("Number of dimensions must be two");
			}
			if ((this->nrow == 0) && (this->ncol == 0)){
				this->nrow = info.shape[0];
				this->ncol = info.shape[1];
				this->AllocateMemory(this->nrow, this->ncol);
			}

			if (info.shape[0] != this->nrow || info.shape[1] != this->ncol) {
				throw std::runtime_error("Input shape does not match matrix dimensions");
			}

			//double* ptr = static_cast<double*>(buf.ptr);

			for (size_t i = 0; i < (size_t)this->nrow; ++i) {
				for (size_t j = 0; j < (size_t)this->ncol; ++j) {
					this->value[i][j] = buf(i,j);// ptr[i * this->ncol + j];
				}
			}
		} /* }}} */
#endif
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
#if HAVE_DEBUG 
			std::cout << "DoubleMatrix: clear memory" << std::endl;
#endif
			DeallocateMemory();
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
		vector<double> operator[](int row){
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

		/* Set dimension */
		void set_dimension(int nrow, int ncol){ /* {{{ */
			this->nrow = nrow;
			this->ncol = ncol;
			AllocateMemory(nrow, ncol);
		} /* }}} */

		/* Set values. */
		void set_value(const vector<vector<double>>& values){ /* {{{ */
			/* Set array's value and matrix size */
			if (this->nrow == 0 && this->ncol == 0){
				this->nrow = (int)values.size();
				this->ncol = (int)values[0].size();

				//if (!this->value){
				if (this->value.size() == 0){
					//cout << "DoubleMatrix::set_value - allocate double** value\n";
					AllocateMemory(this->nrow, this->ncol);
					//this->value = new double*[this->nrow];
					//for (int i=0; i<this->nrow; i++){
					//	this->value[i] = new double[this->ncol];
					//}
				}
				for (int i=0; i<this->nrow; i++){
					for (int j=0; j<this->ncol; j++){
						this->value[i][j] = values[i][j];
					}
				}
			} else{
				if ((int)values.size() !=this->nrow || (int)values[0].size() != this->ncol){
					throw runtime_error("Dimension mismatch!");
				}
				//if (!this->value) {
				if (this->value.size() == 0){
					throw std::runtime_error("Memory for matrix not allocated");
				}
				for (int i=0; i<this->nrow; i++){
					for (int j=0; j<this->ncol; j++){
						this->value[i][j] = values[i][j];
					}
				}
			}
		} /* }}} */
		void set_value(const double A){ /* {{{ */
			int i,j;
			//if (!this->value) {
			if (this->value.size() == 0){
				throw std::runtime_error("Memory for matrix not allocated");
			}

			for (i=0; i<this->nrow; i++){
				for (j=0; j<this->ncol; j++)
					this->value[i][j] = A;
			}
		} /* }}} */
#ifdef HAVE_PYBIND11 
		void set_value_python(py::array_t<double, py::array::c_style>& arr) { /* {{{ */
			py::buffer_info info = arr.request();
			auto buf = arr.mutable_unchecked();

			if (info.ndim != 2) {
				throw std::runtime_error("Number of dimensions must be two");
			}
			if ((this->nrow == 0) && (this->ncol == 0)){
				this->nrow = info.shape[0];
				this->ncol = info.shape[1];
				this->AllocateMemory(this->nrow, this->ncol);
			}

			if (info.shape[0] != this->nrow || info.shape[1] != this->ncol) {
				throw std::runtime_error("Input shape does not match matrix dimensions");
			}

			//double* ptr = static_cast<double*>(buf.ptr);

			for (size_t i = 0; i < (size_t)this->nrow; ++i) {
				for (size_t j = 0; j < (size_t)this->ncol; ++j) {
					this->value[i][j] = buf(i,j);// ptr[i * this->ncol + j];
				}
			}
		} /* }}} */
#endif

		vector<vector<double>> get_value() const{ /* {{{ */
			//if (!this->value) {
			if (this->value.size() == 0){
				throw std::runtime_error("Memory for matrix not allocated");
			}

			vector<vector<double>> result(this->nrow, vector<double>(this->ncol));
			for (int i = 0; i<this->nrow; i++){
				for (int j=0; j<this->ncol; j++){
					result[i][j] = this->value[i][j];
				}
			}
			return result;
		} /* }}} */
#ifdef HAVE_PYBIND11
		py::array_t<double> get_value_python() const{ /* {{{ */
			size_t M = (size_t)this->nrow;
			size_t N = (size_t)this->ncol;

			py::array_t<double, py::array::c_style> arr({M,N});
			auto ra = arr.mutable_unchecked();
			for (size_t i=0; i<M; i++){
				for (size_t j=0; j<N; j++){
					ra(i,j) = this->value[i][j];
				}
			}
			return arr;
		} /* }}} */
		//void get_value_python(py::array_t<double>* arr) { /* {{{ */
		//	size_t M = (size_t)this->nrow;
		//	size_t N = (size_t)this->ncol;

		//	arr = new py::array_t<double>({M,N});
		//	auto ra = arr->mutable_unchecked();
		//	for (size_t i=0; i<M; i++){
		//		for (size_t j=0; j<N; j++){
		//			ra(i,j) = this->value[i][j];
		//		}
		//	}
		//} /* }}} */
#endif

		/* To preventing overloading problem, we try to use "template". */
#if USE_TEMPLATE
		template <typename Template>
		Template get_column(int col) const{ /* {{{ */
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
		} /* }}} */
#else
		double* get_column(int col) const{ /* {{{ */
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
		} /* }}} */
#endif

		double* get_row(int row) const{ /* {{{ */
			if (row > this->nrow || row< 0){
				throw runtime_error("Given row is out of range.");
			}

			double* result = new double[this->ncol];
			for (int i=0; i<this->ncol; i++){
				result[i] = this->value[row][i];
			}

			return result;			
		} /* }}} */
		vector<double> get_row_vector(int row) const{ /* {{{ */
			if (row > this->nrow || row< 0){
				throw runtime_error("Given row is out of range.");
			}

			vector<double> result(this->ncol);
			for (int i=0; i<this->ncol; i++){
				result[i] = this->value[row][i];
			}

			return result;
		} /* }}} */
		
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

// private:
	void AllocateMemory(int nrow, int ncol){ /* {{{ */
		/* Allocate memory*/
		//this->value = new double*[nrow];
		//for (int i=0; i<this->nrow; i++){
		//	this->value[i] = new double[this->ncol];
		//}
		value.resize(nrow, vector<double>(ncol,0.0));
	} /* }}} */

	void DeallocateMemory(){ /* {{{ */
		//int i=0;
		//for (i=0; i<this->nrow; i++){
		//	delete[] value[i];
		//}
		//delete[] value;
		this->value.clear();
	} /* }}} */

	template <class>
    struct always_false : std::false_type {};
};

#endif
