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

```python
import pyseb

% initialize variable
nx = 100 # number of grid
ntime = 365 # number of timestep in day

semic = pyseb.SEMIC()
semic.Initialize(nx)

```

# Check runtime with using OpenMP

![sibma00_openmp_runtime](image/README/sibma00_openmp_runtime.png)

Figure. Runtime of SEMIC depending on number of threads using OpenMP. To operate testing runtime, SIMBA machine is used.
