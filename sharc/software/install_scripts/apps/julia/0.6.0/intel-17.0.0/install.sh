#!/bin/bash
#
# Install Julia 0.6.0 on the ShARC cluster

# Variables
COMPILER=intel
COMPILER_VERS=17.0.0
VERS=0.6.0
MAKEFILE_FPATH=${PWD}/Make.inc
PREFIX="/usr/local/packages/apps/julia/${VERS}/${COMPILER}-${COMPILER_VERS}"

# Error handling
handle_error () {
    errcode=$?
    echo "Error code: $errorcode" 
    echo "Errored command: " echo "$BASH_COMMAND" 
    echo "Error on line: ${BASH_LINENO[0]}"
    exit $errcode
}
trap handle_error ERR

# Load modules
module purge
module load "dev/${COMPILER}-compilers/${COMPILER_VERS}"
module load libs/intel-mkl/2017.0/binary
module load dev/cmake/3.7.1/gcc-4.9.4

# Switch to build directory
pushd "${TMPDIR-/tmp}"

# Clone Julia git repo
if ! [[ -d julia/.git ]]; then
    git clone https://github.com/JuliaLang/julia.git
fi
pushd julia

# Check out a specific version using a tag
git checkout "v${VERS}"

# Modify compilation settings so that use Intel compiler and Intel MKL
cp "${MAKEFILE_FPATH}" .

# Compile
make 
popd

# Delete git history as no longer needed
rm -rf julia/.git

# Copy to installation prefix
mkdir -m 0775 -p "${PREFIX}"
cp -rt julia "${PREFIX}"
