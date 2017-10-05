#!/bin/bash


##########install libtool 2.4.6############
###########################################

# This is a template script for building and installing software on iceberg.
# You should use it to document how you install things.
# You will need to configure any module loads the build needs and then 
# configure the variables for the build.
# This script will then create the directories you need and download and unzip
# the source in to the build dir.

############################# Module Loads ###################################
module load compilers/gcc/4.8.2


############################## Variable Setup ################################
version=2.4.6
prefix=/usr/local/packages6/apps/gcc/4.8.2/flame/0.17.1-openmpi-1.8.8/libtool/$version
build_dir=/data/$USER/FLAME_package/libtool

filename=libtool-$version.tar.gz
baseurl=http://ftpmirror.gnu.org

# Set this to 'sudo' if you want to create the install dir using sudo.
#sudo='sudo'

##############################################################################
# This should not need modifying
##############################################################################

# Create the install directory
if [ ! -d $prefix ]
then
   #$sudo mkdir -p $prefix
   #$sudo chown -R $USER:cs $prefix
   mkdir -p $prefix
   chown -R $USER $prefix
   chgrp -R app-admins $refix
fi 

if [ ! -d $build_dir ]; then
    mkdir -p $build_dir
fi

cd $build_dir

# Download the source
if [ -e $filename ]                                               
then                                                                            
  echo "Install tarball exists. Download not required."                         
else                                                                            
  echo "Downloading source" 
  wget $baseurl/$filename
fi

##############################################################################

##############################################################################
# Installation (Write the install script here)
##############################################################################

tar zxvf $filename

cd libtool-$version

./configure
make
make install prefix=$prefix


##########install libmboard 0.3.1##########
###########################################

# This is a template script for building and installing software on iceberg.
# You should use it to document how you install things.
# You will need to configure any module loads the build needs and then 
# configure the variables for the build.
# This script will then create the directories you need and download and unzip
# the source in to the build dir.

############################# Module Loads ###################################
module load mpi/gcc/openmpi/1.8.8
#make libtool available for libmboard build
export PATH=$PATH://usr/local/packages6/apps/gcc/4.8.2/flame/0.17.1-openmpi-1.8.8/libtool/$version/bin

############################## Variable Setup ################################
version=0.3.1
prefix=/usr/local/packages6/apps/gcc/4.8.2/flame/0.17.1-openmpi-1.8.8/libmboard/$version
build_dir=/data/$USER/FLAME_package

filename=$version.tar.gz
baseurl=https://github.com/FLAME-HPC/libmboard/archive

# Set this to 'sudo' if you want to create the install dir using sudo.
#sudo='sudo'

##############################################################################
# This should not need modifying
##############################################################################

# Create the install directory
if [ ! -d $prefix ]
then
   #$sudo mkdir -p $prefix
   #$sudo chown -R $USER:cs $prefix
   mkdir -p $prefix
   chown -R $USER $prefix
   chgrp -R app-admins $refix
fi 

if [ ! -d $build_dir ]; then
    mkdir -p $build_dir
fi

cd $build_dir

# Download the source
if [ -e $filename ]                                               
then                                                                            
  echo "Install tarball exists. Download not required."                         
else                                                                            
  echo "Downloading source" 
  wget $baseurl/$filename
fi

##############################################################################

##############################################################################
# Installation (Write the install script here)
##############################################################################

tar zxvf $filename

cd $version
#note some additional steps here (possibly)
#cp README.md README.in
#chmod +x create_distribution.should
#./create_distribution.sh


./configure --disable-tests
make
make install prefix=$prefix

##########install xparser 0.17.1##########
###########################################

# This is a template script for building and installing software on iceberg.
# You should use it to document how you install things.
# You will need to configure any module loads the build needs and then 
# configure the variables for the build.
# This script will then create the directories you need and download and unzip
# the source in to the build dir.

############################# Module Loads ###################################
#compilers gcc/4.8.2 & mpi/gcc/openmpi/1.8.8 already loaded


############################## Variable Setup ################################
version=0.17.1
prefix=/usr/local/packages6/apps/gcc/4.8.2/flame/0.17.1-openmpi-1.8.8/xparser/$version
build_dir=/data/$USER/FLAME_package/xparser

filename=$version.tar.gz
baseurl=http://github.com/FLAME-HPC/xparser/archive

# Set this to 'sudo' if you want to create the install dir using sudo.
#sudo='sudo'

##############################################################################
# This should not need modifying
##############################################################################

# Create the install directory
if [ ! -d $prefix ]
then
   #$sudo mkdir -p $prefix
   #$sudo chown -R $USER:cs $prefix
   mkdir -p $prefix
   chown -R $USER $prefix
   chgrp -R app-admins $refix
fi 

if [ ! -d $build_dir ]; then
    mkdir -p $build_dir
fi

cd $build_dir

# Download the source
if [ -e $filename ]                                               
then                                                                            
  echo "Install tarball exists. Download not required."                         
else                                                                            
  echo "Downloading source" 
  wget $baseurl/$filename
fi

##############################################################################

##############################################################################
# Installation (Write the install script here)
##############################################################################

tar zxvf $filename

cd $version

make install prefix=$prefix

bash

