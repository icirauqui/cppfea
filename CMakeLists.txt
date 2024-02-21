cmake_minimum_required(VERSION 3.5)

project(cppfea LANGUAGES CXX)

find_package(MRPT 2.0 COMPONENTS poses gui OPTIONAL_COMPONENTS vision)

# Compile in debug mode
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_VERBOSE_MAKEFILE OFF)


add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/vis)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/fea)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/dataset)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/nlo)


add_executable(cppfea
    src/main.cpp
)

target_link_libraries(cppfea PRIVATE 
    vis
    fea
    dataset
    nlo
    ${MRPT_LIBRARIES}
)