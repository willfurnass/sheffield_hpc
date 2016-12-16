#!/bin/bash

# Install GSL 2.3 on Iceberg

GSL_VERS=2.3
GSL_TARBALL="gsl-${GSL_VERS}.tar.gz"
GSL_TARBALL_URL="http://ftp.heanet.ie/mirrors/gnu/gsl/gsl-2.3.tar.gz"

COMPILER=gcc
COMPILER_VERS=4.9.2

PREFIX="/usr/local/packages6/libs/${COMPILER}/${COMPILER_VERS}/gsl/${GSL_VERS}/"

MODULEFILE="/usr/local/modulefiles/libs/${COMPILER}/${COMPILER_VERS}/gsl/${GSL_VERS}"

BUILD_DIR="${TMPDIR-/tmp}/iceberg/gsl/${GSL_VERS}/${COMPILER}-${COMPILER_VERS}"

WORKERS=${OMP_NUM_THREADS-1}

# Signal handling for failure
handle_error () {
    errcode=$? 
    echo "Error code: $errorcode" 
    echo "Error cmd: " echo "$BASH_COMMAND" 
    echo "Error on line: ${BASH_LINENO[0]}"
    exit $errcode  
}
trap handle_error ERR

# Create and switch to build directory
mkdir -p $BUILD_DIR
pushd $BUILD_DIR

# Download and unpack tarball
[[ -f $GSL_TARBALL ]] || wget $GSL_TARBALL_URL
if ! [[ -f .GSL_TARBALL_unpacked ]]; then
    tar -zxf ${GSL_TARBALL}
    touch .GSL_TARBALL_unpacked 
fi

# Create install and MODULEFILE dirs
for d in $PREFIX $(dirname $MODULEFILE); do
    mkdir -m 2775 -p $d
done

# Activate a compiler
module purge
module load compilers/${COMPILER}/${COMPILER_VERS}

# Build and install
pushd gsl-${GSL_VERS}
./configure --prefix=${PREFIX}
make -j${workers}
make check
make install
popd

# Set permissions and ownership
for d in $PREFIX $(dirname $MODULEFILE); do
    chmod -R g+w $d
    chown -R ${USER}:app-admins $d
done
