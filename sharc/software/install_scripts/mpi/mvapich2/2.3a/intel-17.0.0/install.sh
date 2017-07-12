#!/bin/bash
# Install MVAPICH2 2.3a built with Intel 17.0.0 compilers on the ShARC cluster

##############################################################################
# Error handling
##############################################################################

handle_error () {
    errcode=$? 
    echo "Error: $errorcode" 
    echo "Command: $BASH_COMMAND" 
    echo "Line: ${BASH_LINENO[0]}"
    exit $errcode  
} 
trap handle_error ERR

##############################################################################
# Module loads
##############################################################################

module load dev/intel-compilers/17.0.0
# MVAPICH2 configure script by default errors with "F90 and F90FLAGS are
# replaced by FC and FCFLAGS respectively in this configure, please unset
# F90/F90FLAGS and set FC/FCFLAGS instead" so we now need
for var in F77 F90 F95; do 
    unset $var
done

##############################################################################
# Variable setup
##############################################################################

vers=2.3a
compiler=intel-17.0.0
build_dir="${TMPDIR-/tmp}/${USER}/mvapich2-${vers}"
prefix="/usr/local/packages/mpi/mvapich2/${vers}/${compiler}"
modulefile="/usr/local/modulefiles/mpi/mvapich2/${vers}/${compiler}"
filename="mvapich2-${vers}.tar.gz"
baseurl="http://mvapich.cse.ohio-state.edu/download/mvapich/mv2/"
file_md5="87c3fbf8a755b53806fa9ecb21453445"
workers=16  # for building in parallel

##############################################################################
# Create build, install and modulefile dirs
##############################################################################

[[ -d $build_dir ]] || mkdir -p $build_dir

for d in $prefix $(dirname $modulefile); do 
    mkdir -m 2775 -p $d
    chown -R ${USER}:app-admins $d
done

##############################################################################
# Download source
##############################################################################

cd $build_dir
# Only download if tarball missing or checksum invalid
if ! (echo "${file_md5}  ${filename}" | md5sum -c &> /dev/null); then 
    wget -c ${baseurl}/${filename}
fi

##############################################################################
# Build and install
##############################################################################

if ! [[ -f mvapich2-${vers}/.unpacked ]]; then
    tar -xzf ${filename}
    touch mvapich2-${vers}/.unpacked
fi
cd mvapich2-${vers}

# Enable OmniPath support with --with-device=ch3:psm
# Optional: can specify paths to OmniPath PSM2 headers and libs
# if installed in non-std locations (see --with-psm2-include=path and
# --with-psm2-lib=path 'configure' options)
./configure \
    --prefix=${prefix} \
    --enable-fortran=yes \
    --enable-cxx \
    --with-device=ch3:psm \
    CC=icc CXX=icpc FC=ifort

# SCRAPPY NOTES TO DELETE
#--with-ibverbs-include=??? \
#--with-ibverbs-lib=/usr/lib64/
# FOR OPENMPI: --prefix=${prefix} --with-psm2 CC=icc CXX=icpc FC=ifort

make -j${workers}
make check
make install

###############################################################################
## Download and install examples - DELETE ME
###############################################################################
#
#pushd ${prefix}
#curl -L https://github.com/open-mpi/ompi/archive/v${short_version}.x.tar.gz | tar -zx
#mv ompi-${short_version}.x/examples .
#rm -r ompi-${short_version}.x 
#popd

echo "Next, install modulefile as $modulefile."
