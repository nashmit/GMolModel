cmake_minimum_required(VERSION 3.0)

find_package(Simbody 3.6 CONFIG REQUIRED)
#MESSAGE("Simbody include ")
#MESSAGE(${Simbody_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${Simbody_INCLUDE_DIR})
INCLUDE_DIRECTORIES("/home/nash_mit/git/refactoring/gmolmodel/simbody/SimTKmath/Integrators/src")
LINK_DIRECTORIES(${Simbody_LIBRARY_DIRS})

#MESSAGE("MOLMODEL include ")
#MESSAGE("${Simbody_INCLUDE_DIR/}/..")
INCLUDE_DIRECTORIES(${Simbody_INCLUDE_DIR}/..)
#I should add link_directories for molmodel and add them to TARTGET_LINK!!!

find_package(Eigen3 REQUIRED)
#find_package( PkgConfig )
#pkg_check_modules( EIGEN3 REQUIRED eigen3 )
include_directories( ${EIGEN3_INCLUDE_DIRS} )

#add_definitions(-DWITH_OPENCV2)
add_definitions(-DNUMPY=1)

#change this!!!
#set(mmtk_include "/home/nash_mit/git/refactoring/gmolmodel/mmtk_legacy/Include")
set(mmtk_include "/home/nash_mit/Anaconda/envs/mmtk_environment/include/python2.7/MMTK" )
include_directories( ${mmtk_include} )
MESSAGE(${mmtk_include})
set(mmtk_libs "/home/nash_mit/git/refactoring/gmolmodel/mmtk_legacy/build/lib.linux-x86_64-2.6/MMTK/linux2")
MESSAGE(${mmtk_libs})
link_directories( ${mmtk_libs} )
set(py_include "/home/nash_mit/Anaconda/envs/mmtk_environment/include/python2.7")
include_directories( ${py_include} )

set(numpy_inc "/home/nash_mit/Anaconda/envs/mmtk_environment/lib/python2.7/site-packages/numpy/core/include" )
include_directories( ${numpy_inc} )


#set (Boost_NO_SYSTEM_PATHS ON)
#set (Boost_USE_MULTITHREADED ON)
#set (Boost_USE_STATIC_LIBS ON)
#set (Boost_USE_STATIC_RUNTIME OFF)
#set (BOOST_ALL_DYN_LINK OFF)

find_package(Boost 1.55 REQUIRED) #I'm using 1.58
INCLUDE_DIRECTORIES(${BOOST_INCLUDE_DIR})
LINK_DIRECTORIES(${BOOST_LIBRARY_DIRS})





set(LIBGMOLMODEL_SRC 
	bAddParams.cpp 
	bMainResidue.cpp 
	bMoleculeReader.cpp 
	bSystem.cpp
	)

include_directories(../include/gmolmodel)

add_library(GMOLMODEL_dynamic SHARED ${LIBGMOLMODEL_SRC})
add_library(GMOLMODEL_static  STATIC ${LIBGMOLMODEL_SRC})

# link Simbody libraries
TARGET_LINK_LIBRARIES(GMOLMODEL_dynamic 
		${Simbody_LIBRARIES} 
		${BOOST_LIBRARIES} 
		${Simbody_LIBRARIES} 
		m )
		#${mmtk_libs} )
TARGET_LINK_LIBRARIES(GMOLMODEL_static 
		${Simbody_LIBRARIES} 
		${BOOST_LIBRARIES}
		${Simbody_LIBRARIES} 
		m )
                #${mmtk_libs} )

#set_target_properties(GMOLMODEL_dynamic PROPERTIES OUTPUT_NAME "gmolmodel")
#set_target_properties(GMOLMODEL_dynamic PROPERTIES VERSION 1.2 SOVERSION 1)

#set_target_properties(GMOLMODEL_static PROPERTIES OUTPUT_NAME "gmoldel")

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
