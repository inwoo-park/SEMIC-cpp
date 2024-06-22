# SEMIC-cpp

Wrapping up SEMIC in C++.

# Install CMake

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

# Usage
