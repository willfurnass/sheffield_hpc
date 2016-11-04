#!/bin/bash
#$ -pe openmp 8
#$ -m bea
#$ -M email@sheffield.ac.uk
#$ -j y
#$ -o /home/user/logs/gromacs_2016.1_serial_build.log
#$ -N gromacs_2016_1_serial

################################################
# Install Gromacs 2016.1 (serial) on Iceberg
#
# Will Furnass
# University of Sheffield
################################################

################################################
# Variables
################################################

version=2016.1  
compiler=gcc
compiler_vers=4.9.2  
filename=gromacs-$version.tar.gz
baseurl=ftp://ftp.gromacs.org/pub/gromacs/
build_dir=/scratch/${USER}/gromacs/$version/
inst_dir=/usr/local/packages6/apps/${compiler}/${compiler_vers}/gromacs/${version}-serial
workers=8

# A note re AVX2 (and AVX?) support: the binutils on RHEL6.x is too old to
# support these SIMD flavours in this instance; see
# https://redmine.gromacs.org/issues/1493

################################################
# Signal handling for success and failure
################################################

handle_error () {
    errcode=$? # save the exit code as the first thing done in the trap function 
    echo "error $errorcode" 
    echo "the command executing at the
    time of the error was" echo "$BASH_COMMAND" 
    echo "on line ${BASH_LINENO[0]}"
    # do some error handling, cleanup, logging, notification $BASH_COMMAND
    # contains the command that was being executed at the time of the trap
    # ${BASH_LINENO[0]} contains the line number in the script of that command
    # exit the script or return to try again, etc.
    exit $errcode  # or use some other value or do return instead 
} 
trap handle_error ERR

################################################
# Enable / install dependencies and compilers
################################################

module load compilers/${compiler}/${compiler_vers}
module load compilers/cmake/3.3.0
module load libs/gcc/5.2/boost/1.59

################################################
# Create build and install dirs 
# and set permissions
################################################

mkdir -m 0700 -p ${build_dir}
mkdir -m 2775 -p ${inst_dir}

pushd $build_dir

##################################################
# Download, configure, compile and install Gromacs
##################################################

[[ -d gromacs-$version ]] || curl -L ${baseurl}/${filename} | tar -zx
pushd gromacs-$version
[[ -d build ]] && rm -r build
mkdir build
pushd build
cmake .. \
    -DCMAKE_INSTALL_PREFIX=${inst_dir} \
    -DGMX_FFT_LIBRARY=fftw3 \
    -DGMX_BUILD_OWN_FFTW=ON \
    -DGMX_SIMD=SSE4.1 \
    -DREGRESSIONTEST_DOWNLOAD=ON
    #-DGMX_GPU=ON \
    #-DGMX_MPI=ON \
make -j${workers} 
# Run regression tests
make check
make install

################################################
# Set ownership / permissions
################################################

chown -R ${USER}:app-admins ${inst_dir}
chmod -R g+w ${inst_dir}
