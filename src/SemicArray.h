#ifndef SEMICARRAY_H
#define SEMICARRAY_H

#include <stdexcept>
#include <sstream>
#include <cstring> // for std::memcpy

#ifdef __has_include
#  if __has_include("pybind11/pybind11.h")
#    include "pybind11/pybind11.h"
#    include "pybind11/numpy.h"
#    include "pybind11/stl.h"
#    define HAS_PYBIND11
     namespace py = pybind11;
#  endif
#endif

using namespace std;

class DoubleMatrix {
public:
    DoubleMatrix(size_t rows, size_t cols)
        : nrows(rows), ncols(cols), data(nullptr) {
        allocateMemory();
    }

#ifdef HAS_PYBIND11
    DoubleMatrix(pybind11::array_t<double> array) {
        auto buf = array.request();
        if (buf.ndim != 2) {
            throw std::runtime_error("Number of dimensions must be two");
        }

        nrows = buf.shape[0];
        ncols = buf.shape[1];
        allocateMemory();

        double* ptr = static_cast<double*>(buf.ptr);
        for (size_t i = 0; i < nrows; ++i) {
            std::memcpy(data[i], ptr + i * ncols, ncols * sizeof(double));
        }
    }

    pybind11::array_t<double> toNumpy() const {
        return pybind11::array_t<double>({nrows, ncols}, {ncols * sizeof(double), sizeof(double)}, &data[0][0]);
    }
#endif

    ~DoubleMatrix() {
        freeMemory();
    }

    double& operator()(size_t row, size_t col) {
        if (row >= nrows || col >= ncols) {
            throw std::out_of_range("Index out of range");
        }
        return data[row][col];
    }

    const double& operator()(size_t row, size_t col) const {
        if (row >= nrows || col >= ncols) {
            throw std::out_of_range("Index out of range");
        }
        return data[row][col];
    }

    size_t rows() const { return nrows; }
    size_t cols() const { return ncols; }

    string toString() const {
        ostringstream oss;
        oss << "DoubleMatrix(" << nrows << "x" << ncols << ")";
        return oss.str();
    }

    vector<double> getRow(size_t row) const {
        if (row >= nrows) {
            throw std::out_of_range("Index out of range");
        }
        return std::vector<double>(data[row], data[row] + ncols);
    }

    vector<double> getColumn(size_t col) const {
        if (col >= ncols) {
            throw std::out_of_range("Index out of range");
        }
        std::vector<double> column(nrows);
        for (size_t i = 0; i < nrows; ++i) {
            column[i] = data[i][col];
        }
        return column;
    }

#ifdef HAS_PYBIND11
    void setRow(size_t row, pybind11::array_t<double> array) {
        if (row >= nrows) {
            throw std::out_of_range("Index out of range");
        }
        auto buf = array.request();
        if (buf.ndim != 1 || buf.shape[0] != ncols) {
            throw std::runtime_error("Array must be 1-dimensional and have the correct number of columns");
        }
        double* ptr = static_cast<double*>(buf.ptr);
        std::memcpy(data[row], ptr, ncols * sizeof(double));
    }

    void setColumn(size_t col, pybind11::array_t<double> array) {
        if (col >= ncols) {
            throw std::out_of_range("Index out of range");
        }
        auto buf = array.request();
        if (buf.ndim != 1 || buf.shape[0] != nrows) {
            throw std::runtime_error("Array must be 1-dimensional and have the correct number of rows");
        }
        double* ptr = static_cast<double*>(buf.ptr);
        for (size_t i = 0; i < nrows; ++i) {
            data[i][col] = ptr[i];
        }
    }

    void setMatrix(pybind11::array_t<double> array) {
        auto buf = array.request();
        if (buf.ndim != 2 || buf.shape[0] != nrows || buf.shape[1] != ncols) {
            throw std::runtime_error("Array must be 2-dimensional and match the shape of the matrix");
        }
        double* ptr = static_cast<double*>(buf.ptr);
        for (size_t i = 0; i < nrows; ++i) {
            std::memcpy(data[i], ptr + i * ncols, ncols * sizeof(double));
        }
    }
#endif

private:
    void allocateMemory() {
        data = new double*[nrows];
        for (size_t i = 0; i < nrows; ++i) {
            data[i] = new double[ncols];
        }
    }

    void freeMemory() {
        for (size_t i = 0; i < nrows; ++i) {
            delete[] data[i];
        }
        delete[] data;
    }

    double** data;
    size_t nrows;
    size_t ncols;
};

#endif
