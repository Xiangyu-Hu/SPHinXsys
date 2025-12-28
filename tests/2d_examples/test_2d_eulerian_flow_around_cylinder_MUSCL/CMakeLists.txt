cmake_minimum_required(VERSION 3.10)
project(test_2d_eulerian_flow_around_cylinder_MUSCL LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(test_2d_eulerian_flow_around_cylinder_MUSCL
    2d_eulerian_flow_around_cylinder_MUSCL.cpp
)

target_link_libraries(test_2d_eulerian_flow_around_cylinder_MUSCL
    PRIVATE sphinxsys_2d
)

add_test(NAME test_2d_eulerian_flow_around_cylinder_MUSCL
         COMMAND test_2d_eulerian_flow_around_cylinder_MUSCL)

