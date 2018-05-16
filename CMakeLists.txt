cmake_minimum_required(VERSION 3.10)
project(FilterOpticalFlow)

set(CMAKE_CXX_STANDARD 11)

include_directories(.)

add_executable(FilterOpticalFlow
        src/estimate_velocity.cpp
        include/estimate_velocity/estimate_velocity.h
        src/main.cpp)
