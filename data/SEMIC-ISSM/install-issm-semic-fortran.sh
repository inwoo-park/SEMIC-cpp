#!/bin/bash
: "
Compile original SEMIC with fortran version
"

# Initialize variable
# Before May 26, 2026 ~ Implement SEMIC with C++ in src/c/shared/Semic
COMMIT=893a94277b269f8f5f641281063100a98b8ad8a9
CURRENT_DIR=$(pwd)
ISSM_DIR=${CURRENT_DIR}/ISSM-${COMMIT:0:7}

case "$(hostname)" in
	*gadi*)
		module load python3/3.12.13
		module load openmpi/4.1.3
		module load petsc/3.23

		MPI_DIR="/apps/openmpi/4.1.3"
		;;
	*)
		echo "Error: Not tested at all."
		exit 1;
		;;
esac
PYTHON_DIR="$(dirname $(which python3))/../"
PYTHON_NUMPY_DIR=$(dirname "$(python3 -c 'import numpy; print(numpy.__file__)')")

# Go to ISSM directory
cd ${ISSM_DIR}

# Install SEMIC module
cd ${ISSM_DIR}/externalpackages/semic
./install.sh

# Configure compile
cd ${ISSM_DIR}
if [ ! -f configure ]; then
	autoreconf -ivf
fi
./configure --prefix=${ISSM_DIR} \
	--with-semic-dir=${ISSM_DIR}/externalpackages/semic/install \
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
	LDFLAGS="-L${PETSC_DIR}/lib" CFLAGS="-I${PETSC_DIR}/include"

# Build and install
make -j10
make install
