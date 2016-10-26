cmake_minimum_required(VERSION 2.8)

project(dirichlet_task_solver)

set(CMAKE_CXX_STANDARD 11)

set(CXX_STANDARD_REQUIRED)

SET(CMAKE_CXX_COMPILER mpicxx)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

set(SOURCE_LIB grid.cc gradient_descent.cc mpi_helpers.cc)

set(CXX_COVERAGE_COMPILE_FLAGS "-fopenmp -m64")

set(SOURCE_EXE srcmain.cc)

if(MPI_COMPILE_FLAGS)
  set_target_properties(srcmain PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(srcmain PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

add_library(help_lib STATIC ${SOURCE_LIB} ${MPI_LIBRARIES})

add_executable(srcmain ${SOURCE_EXE})

target_link_libraries(srcmain help_lib)
