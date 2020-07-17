#!/bin/bash

install_dir=/usr/local/packages/libs/hdf5/1.10.4/gcc-8.2.0
mkdir -p $install_dir

tar -xzf hdf5-1.10.4.tar.gz
cd hdf5-1.10.4

module load dev/gcc/8.2
#module load mpi/openmpi/4.0.1/gcc-8.2 

./configure --enable-fortran --enable-cxx --prefix=/usr/local/packages/libs/hdf5/1.10.4/gcc-8.2.0/

#env CC=mpicc FC=mpif90 LDFLAGS=-lnuma CPP=cpp ./configure --enable-fortran --enable-fortran2003 --enable-production --enable-unsupported --enable-parallel --enable-cxx --with-zlib --with-szlib --prefix=$install_dir

make
make check
make install
