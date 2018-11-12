wget https://github.com/lammps/lammps/archive/stable_22Aug2018.tar.gz
tar -xzvf stable_22Aug2018.tar.gz
cd lammps-stable_22Aug2018
mkdir build
cd build
module load dev/cmake/3.7.1/gcc-4.9.4
module load mpi/openmpi/2.0.1/gcc-4.9.4
cmake -C ../cmake/presets/std_nolib.cmake -D CMAKE_INSTALL_PREFIX=/usr/local/packages/apps/lammps/22_08_2018/gcc-4.9.4-openmpi-2.0.1/ BUILD_MPI=yes ../cmake
make
make install
