cmake_minimum_required(VERSION 3.0)

find_package(Simbody 3.6 CONFIG REQUIRED) # you may have to set a path for config file location e.g.  modify the environment variable list specific for find_package..
INCLUDE_DIRECTORIES( ${Simbody_INCLUDE_DIR} )
LINK_DIRECTORIES( ${Simbody_LIB_DIR} )

set( Integrators_PATH 
     "${CMAKE_CURRENT_LIST_DIR}/../simbody/SimTKmath/Integrators/src"
     CACHE FILEPATH "Path to Integrator.h from Simbody..." )
INCLUDE_DIRECTORIES( ${Integrators_PATH} )

set( MOL_MODEL_INCLUDE 
     "${Simbody_INCLUDE_DIR}/.."
     CACHE FILEPATH "header include folder path for Molmodel" ) 
INCLUDE_DIRECTORIES( ${MOL_MODEL_INCLUDE} )

find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIRS} )

add_definitions(-DNUMPY=1)

set(Python_ADDITIONAL_VERSIONS 2)
FIND_PACKAGE(PythonInterp)
FIND_PACKAGE(PythonLibs REQUIRED)
include_directories( ${PYTHON_INCLUDE_DIRS} )
link_directories( ${PYTHON_LIBRARIES} )

set(mmtk_include 
    "${PYTHON_INCLUDE_DIRS}/MMTK" 
    CACHE FILEPATH "Adding header path to MMTK installed location")
include_directories( ${mmtk_include} )

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" 
  OUTPUT_VARIABLE PYTHON_SITE_PATH
  OUTPUT_STRIP_TRAILING_WHITESPACE ) 

set(PYTHON_SITE-PACKAGES_PATH
    ${PYTHON_SITE_PATH} 
    CACHE FILEPATH "side-packages python path")

set(mmtk_libs 
    "${PYTHON_SITE-PACKAGES_PATH}/MMTK/linux2"
    CACHE FILEPATH "Addinng path to MMTK libs folder")
link_directories( ${mmtk_libs} )

execute_process(
  COMMAND ls "${mmtk_libs}"
  OUTPUT_VARIABLE mmtk_LIBRARIES
  OUTPUT_STRIP_TRAILING_WHITESPACE )

set(numpy_include 
    "${PYTHON_SITE-PACKAGES_PATH}/numpy/core/include" 
    CACHE FILEPATH "Header path to numpy installed location")
include_directories( ${numpy_include} )

find_package(Boost 1.55 REQUIRED) #I'm using 1.58
INCLUDE_DIRECTORIES(${BOOST_INCLUDE_DIR})
LINK_DIRECTORIES(${BOOST_LIBRARY_DIRS})

LINK_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/lib)

set(TESTS_SRC 
	TestNative.cpp
	)

include_directories(../include/gmolmodel)

add_executable(TestNative ${TESTS_SRC} )

TARGET_LINK_LIBRARIES( TestNative 
                ${Simbody_LIBRARIES}
                SimTKmolmodel #adding molmodel lib as a dependence.. this needs to be refactor since the name should be changeable
                ${BOOST_LIBRARIES}
                ${PYTHON_LIBRARIES} # do I realy need to link python lib into gmolmodel? it make no sense... but I still don't know verry well code..
                m )

