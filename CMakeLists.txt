cmake_minimum_required(VERSION 3.5)

set(PROJECT_NAME strang-splitting)
project(${PROJECT_NAME})

add_definitions(
	-std=gnu++14
	-Wall
	-pedantic
	)

add_library(strang_lib strang.cpp)

add_executable(simulator simulator.cpp)
target_link_libraries(simulator strang_lib)


############# Eigen ##################
include(FindPkgConfig)
if(PKG_CONFIG_FOUND)
  pkg_search_module(EIGEN REQUIRED eigen3)
  if(EIGEN_FOUND)
    include_directories(${EIGEN_INCLUDE_DIRS})
    message("EIGEN_INCLUDE_DIRS" ${EIGEN_INCLUDE_DIRS})
  endif()
endif()


############# FFTW3 ##################
find_library(FFTW_LIBRARY NAMES fftw3f)
target_link_libraries(strang_lib ${FFTW_LIBRARY})
target_link_libraries(simulator ${FFTW_LIBRARY})


############# Google Test ##################
set(PROJECT_TEST_NAME ${PROJECT_NAME}_test)
find_package(GTest REQUIRED)

# Include gtest library
include_directories(${GTEST_INCLUDE_DIRS})
set(LINK_LIBS ${LINK_LIBS} ${GTEST_BOTH_LIBRARIES})
message("GTEST_INCLUDE_DIRS" ${GTEST_INCLUDE_DIRS})

# somehow this ins necessary...
find_package(Threads)
set(LINK_LIBS ${LINK_LIBS} ${CMAKE_THREAD_LIBS_INIT})


# Unit tests
add_executable(${PROJECT_TEST_NAME} unit_tests.cpp)
#add_dependencies(${PROJECT_TEST_NAME} gtest)

# linking to gtest
target_link_libraries(${PROJECT_TEST_NAME} ${LINK_LIBS})
target_link_libraries(${PROJECT_TEST_NAME} strang_lib)


# make the test runnable with 
# $ make test
enable_testing()
add_test(
	NAME ${PROJECT_TEST_NAME}
	COMMAND ${PROJECT_TEST_NAME}
)