#!/bin/bash
set -eu

HDF5_VERS="1.10.4"
HDF5_SRC_TARBALL="hdf5-${HDF5_VERS}.tar.gz"
HDF5_SRC_TARBALL_MD5="cc451664c06833f9eee5ed29567cd04c"

HDF5_SRC_TARBALL_URL="https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/${HDF5_SRC_TARBALL}"

COMPILER="pgi"
COMPILER_VERS="17.5"

MPI="openmpi"
MPI_VERS="2.0.1"

PREFIX="/usr/local/packages/libs/hdf5/${HDF5_VERS}/${COMPILER}-${COMPILER_VERS}-${MPI}-${MPI_VERS}"
MODULEFILE="/usr/local/modulefiles/libs/hdf5/${HDF5_VERS}/${COMPILER}-${COMPILER_VERS}-${MPI}-${MPI_VERS}"

# Signal handling for failure
handle_error () {
    errcode=$? # save the exit code as the first thing done in the trap function 
    echo "Error: $errcode" 
    echo "Command: $BASH_COMMAND" 
    echo "Line: ${BASH_LINENO[0]}"
    exit $errcode  # or use some other value or do return instead 
}
trap handle_error ERR

# Download and unpack src tarball
[[ -f $HDF5_SRC_TARBALL ]] || wget -L $HDF5_SRC_TARBALL_URL
# MD5 checksums can't be validated - ask about this on HDF5 Forum...
#md5sum ${HDF5_SRC_TARBALL} | grep -q $HDF5_SRC_TARBALL_MD5
if ! [[ -f .hdf5_src_tarball_unpacked ]]; then
    tar -zxf ${HDF5_SRC_TARBALL}
    touch .hdf5_src_tarball_unpacked 
fi

# Create install and modulefile dirs - COMMENTED OUT FOR TESTING
for d in $PREFIX $(dirname $MODULEFILE); do
    mkdir -m 2775 -p $d
done

# Ensure clean environment
module purge

# Load module(s)
module load mpi/openmpi/2.0.1/pgi-17.5
which mpicc 2>&1 > /dev/null

# Build from src and install 
pushd hdf5-${HDF5_VERS}
CFLAGS=-fPIC CC="$(which mpicc)" FC="$(which mpifort)" CPP="$(which cpp)" ./configure --enable-shared --enable-parallel --enable-fortran --enable-fortran2003 --with-zlib --with-szlib --prefix="$PREFIX"
# NB 'CFLAGS=-fPIC' needed - see https://forum.hdfgroup.org/t/link-error-building-parallelhdf5-with-enable-shared/1618
make
make check
make install
popd

## Set permissions and ownership - COMMENTED OUT FOR TESTING
for d in $PREFIX $(dirname $MODULEFILE); do
    chmod -R g+w $d
    chgrp -R hpc_app-admins $d
done
