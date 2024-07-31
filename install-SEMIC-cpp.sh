#!/bin/bash

# get current director
cwd=$(pwd)

if [ ! -d build ]; then
	mkdir build
fi
cd build && cmake ../ -L -DCMAKE_INSTALL_PREFIX=${cwd}/install -DUSE_OPENMP=ON \
	-Dpybind11_DIR=${PYBIND11_DIR} \
	-DUSE_DEBUG=OFF
make && make install
cd ../

#python setup.py build_ext --inplace
