#!/bin/bash

#cd ../molmodel_legacy/simbody
#mkdir build-debug
#cd build-debug

#cmake ..
#make -j4
#sudo make install

#cd ../../
cd ../molmodel_legacy # EU
mkdir build-debug
cd build-debug
cmake ..
make -j12
sudo make install

#cd ../../openmm
#mkdir build-debug
#cd build-debug
#cmake ..
#make -j4
#sudo make install

cd ../../build-debug/


cmake ..
make -j12
