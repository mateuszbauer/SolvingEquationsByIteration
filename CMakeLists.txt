cmake_minimum_required(VERSION 3.18)

project(linear_equations_solving)

set(CMAKE_CXX_STANDARD 17)

add_executable(${PROJECT_NAME} main.cpp matrix.cpp matrix.h)
