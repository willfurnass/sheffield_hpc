#!/bin/bash
#
# The script installs cuDNN 7.6.5.32 for CUDA 10.0 and 9.0 in to the correct ShARC directory.
# The binary can be downloaded from NVIDIA at https://developer.nvidia.com/cudnn
# and need to be saved in $MEDIA_DIR (see below).

set -eu

VERS="7.6.5.32"
VERS_MAJ="${VERS%%.*}"

MEDIA_DIR="/usr/local/media/protected/cuDNN"
PREFIX_BASE="/usr/local/packages/libs/cudnn/$VERS"

for cuda_maj_min_vers in 9.0 10.0 10.1 10.2; do 
    tarball="${MEDIA_DIR}/cudnn-${cuda_maj_min_vers}-linux-x64-v${VERS}.tgz"
    install_dir="${PREFIX_BASE}/binary-cuda-${cuda_maj_min_vers}"
    mkdir -p "$install_dir"
    tar -C "$install_dir" -xvzf "$tarball"

    # Install examples and user manual
    unpack_dir="$(mktemp -d)"
    examples_deb="${MEDIA_DIR}/libcudnn${VERS_MAJ}-doc_${VERS}-1+cuda${cuda_maj_min_vers}_amd64.deb"
    pushd "$unpack_dir"
    ar x "${examples_deb}"
    tar -Jxf data.tar.xz
    mv usr/{src,share} "${install_dir}/cuda"
    popd
done

sudo chgrp -R hpc_app-admins "$PREFIX_BASE"
sudo chmod -R g+w "$PREFIX_BASE"

