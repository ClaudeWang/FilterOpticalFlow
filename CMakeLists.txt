cmake_minimum_required(VERSION 3.10)
project(FilterOpticalFlow)

set(CMAKE_CXX_STANDARD 11)

include_directories(.)

add_executable(FilterOpticalFlow
        EstimateVelocity.cpp
        EstimateVelocity.h
        main.cpp)
