#!/bin/bash
# Install CUDA on the ShARC cluster

# Signal handling for failure
handle_error () {
    errcode=$?
    echo "Error code: $errcode" 
    echo "Errored command: " echo "$BASH_COMMAND" 
    echo "Error on line: ${BASH_LINENO[0]}"
    exit $errcode
}
trap handle_error ERR
set -u

VERSIONS='9.1.85_387.26 8.0.44 7.5.18'
MEDIA_DIR='/usr/local/media/nvidia/'

usage () {
    echo "Usage: $0 [${VERSIONS// /|}]" 1>&2
}

if [ "$#" -ne 1 ]; then
    usage
    exit -1
fi
if [[ ! "$VERSIONS" =~ "$1" ]]; then
    usage
    exit -1
fi

vers="$1"
short_vers="$(echo $1 | cut -d_ -f1)"

pushd $MEDIA_DIR

if [[ $short_vers =~ '9.1' ]]; then
    echo "Checking validity of installers using checksums from https://developer.download.nvidia.com/compute/cuda/9.1/Prod/docs/sidebar/md5sum-3.txt"
    md5sum --check $MEDIA_DIR/cuda-9.1-md5sums.txt
    rc=$?
    if [[ $rc != 0 ]]; then 
        echo 1>&2 'Could not validate MD5 checksums of CUDA 9.1 installer or patches; exiting.'
        exit $rc;
    fi
fi

installer="$MEDIA_DIR/cuda_${vers}_linux.run"
chmod +x "$installer"

prefix="/usr/local/packages/libs/CUDA/${short_vers}/binary"
mkdir -m 2775 -p "$prefix"
chown -R ${USER}:app-admins "$prefix"
chmod -R g+w "$prefix"

echo "Installing CUDA using $installer"
$installer \
    --toolkit --toolkitpath=${prefix}/cuda \
    --samples --samplespath=${prefix}/samples \
    --no-opengl-libs \
    --silent

# Apply patches:
# * cuBLAS Patch Update: This update to CUDA 9.1 includes new GEMM kernels
#   optimized for the Volta architecture and improved heuristics to select GEMM
#   kernels for given input sizes.
# * Patch 2 (Released Feb 27, 2018)     CUDA Compiler Patch Update: This update
#   to CUDA 9.1 includes a bug fix to the PTX assembler (ptxas). The fix resolves
#   an issue when compiling code that performs address calculations using large
#   immediate operands.  
# * Patch 3 (Released Mar 5, 2018)  cuBLAS Patch: This CUDA 9.1 patch includes
#   fixes to GEMM optimizations for convolutional sequence to sequence (seq2seq)
#   models.
if [[ $short_vers =~ '9.1' ]]; then
    for p in 1 2 3; do 
        patch_file="$MEDIA_DIR/cuda_${short_vers}.${p}_linux.run"
        echo "Applying CUDA_related patch $patch_file"
        chmod +x $patch_file
        $patch_file \
           --silent \
           --accept-eula \
           --installdir=${prefix}/cuda
    done 
fi
