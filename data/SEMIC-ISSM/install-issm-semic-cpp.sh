#!/bin/bash
: "
Compile implemented SEMIC in C++ in ISSM
"

# Initialize variable
# Recent version
COMMIT=64da5301fe7c2cb127db5a2fa5899cc31d3560b9
CURRENT_DIR=$(pwd)
ISSM_DIR=${CURRENT_DIR}/ISSM-${COMMIT:1:7}

if [ ! -d ISSM-64da530 ]; then
	git clone https://github.com/ISSMteam/ISSM.git ISSM-${COMMIT:1:7}
fi
cd ${ISSM_DIR}
git checkout ${COMMIT}

case "$(hostname)" in
	*gadi*)
		module load python3/3.12.13
		module load petsc/3.23

		MPI_DIR='/apps/openmpi/4.1.3'
		;;
	*)
		;;
esac
PYTHON_DIR="$(dirname $(which python3))/../"
PYTHON_NUMPY_DIR=$(dirname "$(python3 -c 'import numpy; print(numpy.__file__)')")

# Configure compile
cd ${ISSM_DIR}
if [ ! -f configure ]; then
	autoreconf -ivf
fi
./configure --prefix=${ISSM_DIR} \
	--with-semic=${ISSM_DIR}/externalpackages/semic/install \
	--with-python-dir="${PYTHON_DIR}" \
	--with-python-numpy-dir="${PYTHON_NUMPY_DIR}" \
	--with-mpi-include="${MPI_DIR}/include" \
	--with-mpi-libflags="-L${MPI_DIR}/lib -lmpi " \
	--with-metis-dir="${PETSC_DIR}" \
	--with-parmetis-dir="${PETSC_DIR}" \
	--with-blas-lapack-dir="${PETSC_DIR}" \
	--with-scalapack-dir="${PETSC_DIR}" \
	--with-mumps-dir="${PETSC_DIR}" \
	--with-petsc-dir="${PETSC_DIR}" \
	LDFLAGS="-L${PETSC_DIR}/lib"

# Build and install
make -j10
make install
