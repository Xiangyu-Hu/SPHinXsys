cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(test_shock_tube LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_NAME "test_shock_tube")

AUX_SOURCE_DIRECTORY(. DIR_SRCS)

ADD_EXECUTABLE(${PROJECT_NAME} ${DIR_SRCS})

target_link_libraries(${PROJECT_NAME} PRIVATE sphinxsys_2d)

find_package(GTest REQUIRED)

target_link_libraries(${PROJECT_NAME} PRIVATE GTest::gtest GTest::gtest_main)

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
