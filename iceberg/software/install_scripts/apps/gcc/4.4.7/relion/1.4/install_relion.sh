#!/bin/bash -e

module purge

version=1.4
install_dir=/usr/local/packages6/apps/gcc/4.4.7/relion/$version

rm -rf /scratch/te1st/relion
mkdir -p /scratch/te1st/relion
cd /scratch/te1st/relion

wget http://www2.mrc-lmb.cam.ac.uk/groups/scheres/1sep15/relion-1.4.tar.bz2
tar -xjf relion-1.4.tar.bz2

cd relion-1.4

module load mpi/gcc/openmpi/1.8.8
./configure --enable-mpi --prefix=$install_dir
make
make install


