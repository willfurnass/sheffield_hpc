#!/bin/bash
#$ -N openblas_iceberg
#$ -j y
#$ -pe openmp 4
export OMP_NUM_THREADS=4

OPENBLAS_VERS=0.2.19
OPENBLAS_TARBALL="v${OPENBLAS_VERS}.tar.gz"
OPENBLAS_TARBALL_URL="http://github.com/xianyi/OpenBLAS/archive/${OPENBLAS_TARBALL}"

COMPILER=gcc
COMPILER_VERS=4.9.2

BUILD_DIR="${TMP_DIR-/tmp}/openblas/${OPENBLAS_VERS}/${COMPILER}-${COMPILER_VERS}"

PREFIX="/usr/local/packages6/libs/${COMPILER}/${COMPILER_VERS}/openblas/${OPENBLAS_VERS}/"

MODULEFILE="/usr/local/modulefiles/libs/${COMPILER}/${COMPILER_VERS}/openblas/${OPENBLAS_VERS}"

# Signal handling for failure
handle_error () {
    errcode=$? 
    echo "Error code: $errorcode" 
    echo "Error cmd: " echo "$BASH_COMMAND" 
    echo "Error on line: ${BASH_LINENO[0]}"
    exit $errcode  
}
trap handle_error ERR

echo "Installing on: ${HOSTNAME} ($(grep -m1  Xeon /proc/cpuinfo | cut -d: -f2))"

# Make and switch to build dir
mkdir -p $BUILD_DIR
pushd $BUILD_DIR

# Download and unpack tarball
[[ -f $OPENBLAS_TARBALL ]] || wget $OPENBLAS_TARBALL_URL
if ! [[ -f .OPENBLAS_TARBALL_unpacked ]]; then
    tar -zxf ${OPENBLAS_TARBALL}
    touch .OPENBLAS_TARBALL_unpacked 
fi

# Create install and modulefile dirs
for d in $PREFIX $(dirname $MODULEFILE); do
    mkdir -m 2775 -p $d
done

# Activate a compiler
module purge
module load compilers/${COMPILER}/${COMPILER_VERS}

# Build and install
pushd OpenBLAS-${OPENBLAS_VERS}
make -j $OMP_NUM_THREADS
make PREFIX=${PREFIX} install
popd

# Set permissions and ownership
for d in $PREFIX $(dirname $MODULEFILE); do
    chmod -R g+w $d
    chown -R ${USER}:app-admins $d
done
