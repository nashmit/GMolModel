#!/bin/bash

sudo apt-get install gfortran
sudo apt-get install libglfw3-dev
sudo apt-get install freeglut3-dev
sudo apt-get install libglew-dev
sudo apt-get install libxmu-dev
sudo apt-get install libeigen3-dev
sudo apt-get install doxygen
sudo apt-get install subversion
sudo apt-get install libblas-dev liblapack-dev
sudo apt-get install libboost-all-dev
sudo apt-get install swig
sudo apt-get install ocl-icd-opencl-dev

cd ../molmodel_legacy/simbody
mkdir build-debug
cd build-debug

cmake ..
make -j4
sudo make install

cd ../../
mkdir build-debug
cd build-debug
cmake ..
make -j4
sudo make install

cd ../../openmm
mkdir build-debug
cd build-debug
cmake ..
make -j4
sudo make install

cd ../../build-debug/


cmake .. -DMakeDownload_mmtk_environment=true -DInstallEnvironment=true -DAnacondaEnvironmentPath="/home/laurentiu/anaconda2/envs/"

