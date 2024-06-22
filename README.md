# SEMIC-cpp

Wrapping up SEMIC in C++.

# Install

## Required Packages in python

* pybind11
* numpy
* matplotlib

## Install

* To install the `pySEMIC.cpp` for running semic in cpp, use following command.

* This work is only tested in `linux` machine, `not window` machine.

```bash
cwd = $(pwd)
mkdir build
cd build
cmake ../ -L -DCMAKE_INSTALL_PREFIX=${cwd}/../install
make 
make install
```

## Window system

* Require Miscrotsoft Visual C++ 14.0 or greater.
* Download https://visualstudio.microsoft.com/visual-cpp-build-tools/
* Install `developement` verison of current project as follows:
* Window 10 > Visual Studio 2019
* Window 11 > Visual Studio 2019 or 2022

```
pip install -e .
```

# Usage
