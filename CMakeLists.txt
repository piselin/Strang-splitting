cmake_minimum_required(VERSION 3.5)

set(PROJECT_NAME strang-splitting)
project(${PROJECT_NAME})

set(target_name simulator)

add_definitions(
	-std=gnu++14
	-Wall
	-O3
)

# Add this definition if the FFTW3 backend should be used
# without this, the default backend kissFFT will be used.
add_definitions(-DEIGEN_FFTW_DEFAULT)

add_executable(${target_name} simulator.cpp)


############# Eigen ##################
find_package(Eigen3 REQUIRED QUIET)
include(FindPkgConfig)
if(PKG_CONFIG_FOUND)
	pkg_search_module(EIGEN REQUIRED eigen3)
	if(EIGEN_FOUND)
		include_directories(${EIGEN_INCLUDE_DIRS})
		message("EIGEN_INCLUDE_DIRS " ${EIGEN_INCLUDE_DIRS})
	endif()
endif()

############# FFTW3 ##################

# the default should be "OFF"
option(FFTW_SINGLE_PRECISION "Use FFTW single precision" OFF)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(FFTW REQUIRED)
if(FFTW_FOUND)
	message("FFTW3_INCLUDE_DIRS " ${FFTW_INCLUDES})

	if(FFTW_SINGLE_PRECISION)
		message("Using single point precission " ${FFTW_LIB_FLOAT})
		set(CHOSEN_FFTW_LIB ${FFTW_LIB_FLOAT})
	else()
		message("Using double point precission " ${FFTW_LIB})
		set(CHOSEN_FFTW_LIB ${FFTW_LIB})
	endif()

	include_directories(${FFTW_INCLUDES_DIRS})
	target_link_libraries(${target_name} ${CHOSEN_FFTW_LIB})

endif()

############# HDF5 ##################
find_package(HDF5 COMPONENTS CXX)
if (HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIR})
  message("HDF5_INCLUDE_DIRS " ${HDF5_INCLUDE_DIR})
  set(LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})
endif()


############# Google Test ##################
set(PROJECT_TEST_NAME ${PROJECT_NAME}-test)
find_package(GTest REQUIRED)

# Include gtest library
include_directories(${GTEST_INCLUDE_DIRS})
set(LINK_LIBS ${LINK_LIBS} ${GTEST_BOTH_LIBRARIES})
message("GTEST_INCLUDE_DIRS " ${GTEST_INCLUDE_DIRS})

# Necessary for Google Test
find_package(Threads)
set(LINK_LIBS ${LINK_LIBS} ${CMAKE_THREAD_LIBS_INIT})

# Unit tests
add_executable(${PROJECT_TEST_NAME} unit_tests.cpp)

# linking to gtest
target_link_libraries(${PROJECT_TEST_NAME} ${LINK_LIBS})
target_link_libraries(${PROJECT_TEST_NAME} ${CHOSEN_FFTW_LIB})

# make the test runnable with 
# $ make test
enable_testing()
add_test(
	NAME ${PROJECT_TEST_NAME}
	COMMAND ${PROJECT_TEST_NAME}
)


############## Link everything together ##################
target_link_libraries(${target_name} ${LINK_LIBS})
