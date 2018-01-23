#!/bin/bash

petsc_vers=3.8.3
petsc_tarball="petsc-${petsc_vers}.tar.gz"
petsc_tarball_url="http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/${petsc_tarball}"

compiler_info='intel-17.0.0-openmpi-2.0.1'
compiler_modfile='mpi/openmpi/2.0.1/intel-17.0.0'
blas_modfile='libs/intel-mkl/2017.0/binary'

prefix="/usr/local/packages/libs/petsc/${petsc_vers}/${compiler_info}"

modulefile="/usr/local/modulefiles/libs/petsc/${petsc_vers}/${compiler_info}"

# Need to build in a dir on a shared filesystem otherwise the MPI tests won't work
unpack_dir="/data/$USER/building/${SGE_CLUSTER_NAME-$(localhost)}_petsc/${petsc_vers}"
build_dir="${unpack_dir}/${compiler_info}"

# Signal handling for failure
handle_error () {
    errcode=$? # save the exit code as the first thing done in the trap function 
    echo "Error: $errorcode" 
    echo "Command: $BASH_COMMAND" 
    echo "Line: ${BASH_LINENO[0]}"
    exit $errcode  # or use some other value or do return instead 
}
trap handle_error ERR

# Download and unpack tarball
mkdir -p $unpack_dir
pushd $unpack_dir
wget -N $petsc_tarball_url
if ! [[ -f .petsc_tarball_unpacked ]]; then
    tar -zxf ${petsc_tarball}
    touch .petsc_tarball_unpacked 
fi

# Create install and modulefile dirs
for d in $prefix $(dirname $modulefile); do
    #mkdir -m 2775 -p $d
    mkdir -p $d
done

# Activate a (MPI) compiler and BLAS
module purge
module load ${blas_modfile}
module load ${compiler_modfile}

# Configure, build and install 
mkdir -p $build_dir
cp -Tra petsc-${petsc_vers} ${build_dir}
pushd "${build_dir}"

# _Should_ automatically find the active MPI 
python2 ./configure \
    --prefix=${prefix} \
    --with-debugging=0 \
    --with-blaslapack-dir=$MKLROOT \
    --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpifort \
    CFLAGS=$CFLAGS \
    CPPFLAGS=$CPPFLAGS \
    LDFLAGS=$LDFLAGS \
    COPTFLAGS='-O3 -march=native -mtune=native' CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native'
make clean
make MAKE_NP=${NSLOTS-1} all 
make MAKE_NP=${NSLOTS-1} test 
make install
popd

## Set permissions and ownership
#for d in $prefix $(dirname $modulefile); do
#    chmod -R g+w $d
#    chown -R ${USER}:app-admins $d
#done
