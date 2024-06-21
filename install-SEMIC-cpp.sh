#!/bin/bash

# get current director
cwd=$(pwd)

rm -rf build
mkdir build
cd build && cmake ../ -L -DCMAKE_INSTALL_PREFIX=${cwd}/install -DUSE_OPENMP=OFF
make && make install
