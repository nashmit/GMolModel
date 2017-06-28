#!/bin/bash

sudo apt-get install gfortran
sudo apt-get install libglfw3-dev
sudo apt-get install freeglut3-dev
sudo apt-get install freeglut3-dev
sudo apt-get install libxmu-dev
sudo apt-get install libeigen3-dev
sudo apt-get install doxygen
sudo apt-get install libblas-dev liblapack-dev
sudo apt-get install libboost-all-dev

cmake .. -DMakeDownload_mmtk_environment=true -DInstallEnvironment=true -DAnacondaEnvironmentPath="/home/nash_mit/Anaconda/envs/"

