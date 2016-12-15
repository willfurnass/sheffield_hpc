#!/bin/bash
# Install CASTEP 16.11 on the Iceberg cluster

##############################################################################
# Error handling
##############################################################################

handle_error () {
    errcode=$? # save the exit code as the first thing done in the trap function 
    echo "Error code: $errorcode" 
    echo "Errored command: " echo "$BASH_COMMAND" 
    echo "Error line: ${BASH_LINENO[0]}"
    exit $errcode  
} 
trap handle_error ERR

##############################################################################
# Variables
##############################################################################

vers="16.11"
compiler="intel"
compiler_short_vers="15"
compiler_vers="15.0.3"
tarball_path="/usr/local/media/protected/CASTEP/${vers}/CASTEP-${vers}.tar.gz"
base_prefix="/usr/local/packages6/apps/"
prefix="${base_prefix}/${compiler}/${compiler_short_vers}/castep/${vers}"
serial_build_dir="/data/$USER/iceberg/castep/${vers}/${compiler}-${compiler_short_vers}/serial"
parallel_build_dir="/data/$USER/iceberg/castep/${vers}/${compiler}-${compiler_short_vers}/parallel"
#export OMP_NUM_THREADS=4
#workers=$OMP_NUM_THREADS
workers=1

##############################################################################
# Create dirs
##############################################################################

mkdir -m 0700 -p $serial_build_dir $parallel_build_dir
mkdir -m 2750 -p $prefix

##############################################################################
# Build and install serial version
##############################################################################

module purge
module load compilers/${compiler}/${compiler_vers}
module load libs/binlibs/intel-mkl/11.2.3

tar -xzf ${tarball_path} -C ${serial_build_dir}
pushd ${serial_build_dir}/CASTEP-${vers}

make clean
make -j $workers INSTALL_DIR=$prefix FFT=mkl MATHLIBS=mkl10
make -j $workers INSTALL_DIR=$prefix FFT=mkl MATHLIBS=mkl10 install install-tools

##############################################################################
# Build and install parallel version
##############################################################################

module load mpi/intel/openmpi/1.10.0

tar -xzf ${tarball_path} -C ${parallel_build_dir}
pushd ${parallel_build_dir}/CASTEP-${vers}

if [[ $compiler == "intel" ]] && [[ $compiler_vers =~ ^15\. ]]; then
    # Workaround for bug described at http://www.cmth.ph.ic.ac.uk/computing/software/castep.html
    sed 's/-static-intel/-shared-intel/' obj/platforms/linux_x86_64_ifort15.mk -i
fi

make clean
make -j $workers COMMS_ARCH=mpi  FFT=mkl MATHLIBS=mkl10
mv ./obj/linux_x86_64_*/castep.mpi $prefix

##############################################################################
# Set permissions
##############################################################################

chmod -R g-w $prefix
chmod -R o-rwx $prefix
chgrp -R castep $prefix
