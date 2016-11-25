#!/bin/bash

#The script installs cuDNN 5.1 for CUDA 7.5 and CUDA 8.0 in to the correct iceberg directory

#The binaries can be downloaded from Nvidia at https://developer.nvidia.com/cudnn

#Binaries on Iceberg are located at /usr/local/media/protected/cuDNN

#We assume that the .tgz fils are in current working directory

#Install cuDNN 5.0 for CUDA 7.5
install_dir=/usr/local/packages6/libs/binlibs/cudnn/5.1-cuda-7.5
mkdir -p $install_dir
tar -C $install_dir -xvzf ./cudnn-7.5-linux-x64-v5.1-ga.tgz 

#Install cuDNN 5.0 for CUDA 8.0
install_dir=/usr/local/packages6/libs/binlibs/cudnn/5.1-cuda-8.0
mkdir -p $install_dir
tar -C $install_dir -xvzf ./cudnn-8.0-linux-x64-v5.1-ga.tgz 

