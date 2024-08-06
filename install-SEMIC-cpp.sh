#!/bin/bash

# get current director
cwd=$(pwd)

if [ ! -d build ]; then
	mkdir build
fi
cd build && cmake ../ \
	--fresh \
	-L -DCMAKE_INSTALL_PREFIX=${cwd}/install -DUSE_OPENMP=ON \
	-DUSE_DEBUG=OFF
make && make install
cd ../

#-Dpybind11_DIR=${PYBIND11_DIR} \
#python setup.py build_ext --inplace
