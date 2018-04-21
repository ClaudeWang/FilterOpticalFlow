cmake_minimum_required(VERSION 3.5.1)
project(FilterOpticalFlow)

set(CMAKE_CXX_STANDARD 11)

add_executable(FilterOpticalFlow main.cpp EstimateVelocity.cpp EstimateVelocity.h)