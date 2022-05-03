cmake_minimum_required(VERSION 3.18)

project(test)

set(CMAKE_CXX_STANDARD 17)

add_executable(test main.cpp matrix.cpp matrix.h)
