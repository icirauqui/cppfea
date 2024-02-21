cmake_minimum_required(VERSION 3.5)

project(nlo LANGUAGES CXX)

find_package(MRPT 2.0 COMPONENTS poses gui OPTIONAL_COMPONENTS vision)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_VERBOSE_MAKEFILE OFF)


add_library(nlo
    levenberg_marquardt.cpp
    levenberg_marquardt.hpp
)

target_link_libraries(nlo PRIVATE 
    fea
    ${MRPT_LIBRARIES}
)













