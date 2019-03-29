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

VERSIONS='10.0.130_410.48 9.1.85_387.26 9.0.176_384.81 8.0.44 7.5.18'
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
# The short version is everything up to but not including the first underscore
short_vers="${vers%%_*}"

pushd $MEDIA_DIR

if [[ $short_vers =~ '^9\.' ]]; then
    echo "Checking validity of installers using checksums from NVIDIA site"
    md5sum --check $MEDIA_DIR/cuda-${short_vers}-md5sums.txt
    rc=$?
    if [[ $rc != 0 ]]; then 
        echo 1>&2 "Could not validate MD5 checksums of CUDA ${short_vers} installer or patches; exiting."
        exit $rc;
    fi
fi

installer="$MEDIA_DIR/cuda_${vers}_linux.run"
echo $short_vers
if [[ $short_vers =~ ^9.0[^0-9] ]]; then
    installer="$MEDIA_DIR/cuda_${vers}_linux-run"
fi
chmod +x "$installer"

prefix="/usr/local/packages/libs/CUDA/${short_vers}/binary"
mkdir -m 2775 -p "$prefix"
chown -R ${USER}:hpc_app-admins "$prefix"
chmod -R g+w "$prefix"

echo "Installing CUDA using $installer"
$installer \
    --toolkit --toolkitpath=${prefix}/cuda \
    --samples --samplespath=${prefix}/samples \
    --no-opengl-libs \
    --silent

# Apply patches:
if [[ $short_vers =~ ^9.1[^0-9] ]]; then
    for p in 1 2 3; do 
        patch_file="$MEDIA_DIR/cuda_${short_vers}.${p}_linux.run"
        echo "Applying CUDA_related patch $patch_file"
        chmod +x $patch_file
        $patch_file \
           --silent \
           --accept-eula \
           --installdir=${prefix}/cuda
    done 
elif [[ $short_vers =~ ^9.0[^0-9] ]]; then
    for p in 1 2 3; do 
        patch_file="$MEDIA_DIR/cuda_${short_vers}.${p}_linux-run"
        echo "Applying CUDA_related patch $patch_file"
        chmod +x $patch_file
        $patch_file \
           --silent \
           --accept-eula \
           --installdir=${prefix}/cuda
    done 
fi
