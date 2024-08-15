#!/bin/bash

CC=gcc
CXX=g++
if [ $(hostname) == 'simba20' ] || [ $(hostname) == 'simba00' ]; then
	CC=icc
	CXX=icpc
fi

# get current director
cwd=$(pwd)

if [ ! -d build ]; then
	mkdir build
fi

options=''
if [ -f build/CMakeLists.txt ]; then
	options="$options --fresh"
fi

echo $options

cd build && cmake ../ $options \
	-L -DCMAKE_INSTALL_PREFIX=${cwd}/install -DUSE_OPENMP=ON \
	-DUSE_DEBUG=OFF \
	-DCMAKE_C_COMPILER=${CC} \
	-DCMAKE_CXX_COMPILER=${CXX}

make && make install
cd ../

#-Dpybind11_DIR=${PYBIND11_DIR} \
#python setup.py build_ext --inplace
