#!/bin/bash
# Install CUDA on the ShARC cluster

# Signal handling for failure
handle_error () {
    errcode=$?
    echo "Error code: $errorcode" 
    echo "Errored command: " echo "$BASH_COMMAND" 
    echo "Error on line: ${BASH_LINENO[0]}"
    exit $errcode
}
trap handle_error ERR

VERSIONS="9.1.85_387.26 8.0.44 7.5.18"

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

version="$1"
short_version="$(echo $1 | cut -d_ -f1)"

installer="/usr/local/media/nvidia/cuda_${version}_linux.run"
chmod +x "$installer"

prefix="/usr/local/packages/libs/CUDA/${short_version}/binary"
mkdir -m 2775 -p "$prefix"
chown -R ${USER}:app-admins "$prefix"
chmod -R g+w "$prefix"

$installer \
    --toolkit --toolkitpath=${prefix}/cuda \
    --samples --samplespath=${prefix}/samples \
    --no-opengl-libs \
    --silent
