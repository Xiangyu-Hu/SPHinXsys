cmake_minimum_required(VERSION 3.10)
project(test_2d_energy_gradient LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GTest REQUIRED)

add_executable(test_2d_energy_gradient
    2d_energy_gradient.cpp
)

target_link_libraries(test_2d_energy_gradient
    PRIVATE sphinxsys_2d
    PRIVATE GTest::gtest_main
    PRIVATE GTest::gtest
)

add_test(NAME test_2d_energy_gradient COMMAND test_2d_energy_gradient)

