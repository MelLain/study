cmake_minimum_required(VERSION 2.8)

project(dirichlet_task_solver)

set(CMAKE_CXX_COMPILER "g++")

set(SOURCE_LIB grid.cc gradient_descent.cc)

set(CMAKE_CXX_STANDARD 11)

set(CXX_STANDARD_REQUIRED)

set(CXX_COVERAGE_COMPILE_FLAGS "-fopenmp -m64")

set(SOURCE_EXE srcmain.cc)

add_library(help_lib STATIC ${SOURCE_LIB})

add_executable(srcmain ${SOURCE_EXE})

target_link_libraries(srcmain help_lib)
