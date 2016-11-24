#!/bin/bash
################################################################################
# Install GULP 4.4 (serial and MPI builds) on the ShARC cluster
################################################################################

################################################################################
# Variables
################################################################################
vers=4.4
tarball_path=/usr/local/media/protected/GULP/${vers}/gulp-${vers}.tgz
compiler_vers=17.0.0
openmpi_vers=2.0.1

serial_build_dir=${TMPDIR-/tmp}/$USER/gulp/$vers/intel-${compiler_vers}/serial
serial_build_config=$(realpath getmachine_serial)
serial_prefix=/usr/local/packages/apps/gulp/${vers}/intel-${compiler_vers}

mpi_build_dir=${TMPDIR-/tmp}/$USER/gulp/$vers/intel-${compiler_vers}/mpi
mpi_build_config=$(realpath getmachine_mpi)
mpi_prefix=/usr/local/packages/apps/gulp/${vers}/intel-${compiler_vers}-openmpi-${openmpi_vers}

workers=16

################################################################################
# Error handling
################################################################################
handle_error () {
    errcode=$? 
    echo "Error code: $errorcode" 
    echo "Errored command: " echo "$BASH_COMMAND" 
    echo "Error on line: ${BASH_LINENO[0]}"
    exit $errcode  
}
trap handle_error ERR

################################################################################
# Serial build
################################################################################
module purge
module load dev/intel-compilers/17.0.0
module load libs/intel-mkl/2017.0/binary

mkdir -m 0700 -p $serial_build_dir
mkdir -m 2770 -p $serial_prefix

tar -zxf $tarball_path -C $serial_build_dir
pushd $serial_build_dir/gulp-${vers}/Src
# Create backup of config file
[[ -f getmachine.orig ]] || cp getmachine{,.orig}
# Copy config file for serial build in place
cp $serial_build_config getmachine
make -j $workers

cp -r $serial_build_dir/gulp-${vers}/* $serial_prefix/
mkdir $serial_prefix/bin
ln -sr $serial_prefix/Src/gulp -t $serial_prefix/bin/

################################################################################
# MPI build
################################################################################
module load mpi/openmpi/2.0.1/intel-17.0.0

mkdir -m 0700 -p $mpi_build_dir
mkdir -m 2770 -p $mpi_prefix

tar -zxf $tarball_path -C $mpi_build_dir
pushd $mpi_build_dir/gulp-${vers}/Src
# Create backup of config file
[[ -f getmachine.orig ]] || cp getmachine{,.orig}
# Copy config file for MPI build in place
cp $mpi_build_config getmachine
make -j $workers

cp -r $mpi_build_dir/gulp-${vers}/* $mpi_prefix/
mkdir $mpi_prefix/bin
ln -sr $mpi_prefix/Src/gulp -t $mpi_prefix/bin/

################################################################################
# Restrict permissions
################################################################################
chmod -R o-rwx $serial_prefix $mpi_prefix
echo "Now, run the following: sudo chgrp -R gulp $serial_prefix $mpi_prefix"
