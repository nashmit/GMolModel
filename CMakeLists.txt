cmake_minimum_required(VERSION 3.0)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

project(gmolmodel VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -fPIC -g")

set( MakeDownload_mmtk_environment 
     FALSE
     CACHE BOOL "Should I download the anaconda's mmtk_environment?" )

if( MakeDownload_mmtk_environment )
        set( wget_parameter "ftp://mybioinformaticslab.xyz/mmtk_environment.tar.gz" )
        execute_process(
                COMMAND wget ${wget_parameter}
                OUTPUT_VARIABLE statusWGET )
	MESSAGE( STATUS  "wget status: ${wget_parameter}" )
	execute_process(
		COMMAND tar -xzvf "mmtk_environment.tar.gz"
		OUTPUT_VARIABLE statusUnzip )
	MESSAGE( STATUS "tar status ${statusUnzip}" )
endif( MakeDownload_mmtk_environment )

set( InstallEnvironment 
     false 
     CACHE BOOL "Install mmtk anaconda environment" )
set( AnacondaEnvironmentPath 
     ""
     CACHE FILEPATH "envs path for Anaconda: e.g. :/home/nash_mit/Anaconda/envs/" )
    
if( InstallEnvironment AND AnacondaEnvironmentPath ) 
	execute_process(
		COMMAND cp -RfHp "mmtk_environment" ${AnacondaEnvironmentPath}
		OUTPUT_VARIABLE statusCP )
	MESSAGE( STATUS "Done with installing mmtk_environment ${statusCP} after this you should: source activate mmtk_environment" )
elseif( InstallEnvironment AND  NOT AnacondaEnvironmentPath )
	MESSAGE( FATAL_ERROR "If you want to install Anaconda environment you must set AnacondaEnvironmentPath " )
endif( InstallEnvironment AND AnacondaEnvironmentPath )

if( NOT InstallEnvironment )
	add_subdirectory(src)
endif( NOT InstallEnvironment )

