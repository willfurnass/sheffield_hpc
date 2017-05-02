#!/bin/bash

# This is a template script for building and installing software on iceberg.
# You should use it to document how you install things.
# You will need to configure any module loads the build needs and then 
# configure the variables for the build.
# This script will then create the directories you need and download and unzip
# the source in to the build dir.

############################# Module Loads ###################################
module load dev/cmake/3.7.1/gcc-4.9.4

############################## Variable Setup ################################
version=2.1.0
prefix=/usr/local/packages/apps/ANTs/$version
package_dir=/scratch/$USER/ANTs_package
source_dir=$package_dir/ANTs
build_dir=$package_dir/ANTs_build

filename=ANTs.git
baseurl=https://github.com/stnava

# Set this to 'sudo' if you want to create the install dir using sudo.
#sudo='sudo'

##############################################################################
# This should not need modifying
##############################################################################

# Create the install directory
if [ ! -d $prefix ]
then
   mkdir -p $prefix
   chown -R $USER $prefix
   chgrp -R app-admins $prefix
fi 

#create package directory
if [ ! -d $package_dir ]; then
    mkdir -p $package_dir
fi

#create source directory
#if [ ! -d $source_dir ]; then
#    mkdir -p $source_dir
#fi

#create build directory
if [ ! -d $build_dir ]; then
    mkdir -p $build_dir
fi

cd $package_dir

# Download the source
if [ -e $filename ]                                               
then                                                                            
  echo "Install tarball exists. Download not required."                         
else                                                                            
  echo "Downloading source" 
  git clone $baseurl/$filename
fi

##############################################################################

##############################################################################
# Installation (Write the install script here)
##############################################################################

#configure uild directory using cmake
cmake -B$build_dir -H$source_dir

cd $build_dir

make

#install build on HPC


#mv ./STAR /usr/local/packages6/apps/gcc/5.2/STAR/2.5.0c/
cp -R $build_dir/* $prefix

bash

