cmake_minimum_required(VERSION 3.5)

project(dataset LANGUAGES CXX)

set(CMAKE_VERBOSE_MAKEFILE OFF)


add_library(dataset
    dataset.cpp
    dataset.hpp
)