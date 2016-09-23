cmake_minimum_required(VERSION 2.8)

SET(GCC_COVERAGE_COMPILE_FLAGS "-std=c++11 -stdlib=libc++")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

project(lab_algorithmic_one)

set(SOURCE_EXE srcmain.cc)

set(SOURCE_LIB helper.cc)

add_library(help_lib STATIC ${SOURCE_LIB})

add_executable(srcmain ${SOURCE_EXE})

target_link_libraries(srcmain help_lib)
