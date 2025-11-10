cmake_minimum_required(VERSION 3.15)
project(Chislaki CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

add_subdirectory(src)
add_subdirectory(Quadratic_Function_Minimization_using_Gradient_and_Coordinate_Descent_Methods)
