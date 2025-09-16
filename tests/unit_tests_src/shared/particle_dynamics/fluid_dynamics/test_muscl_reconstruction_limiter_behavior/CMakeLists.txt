cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(test_muscl_reconstruction_limiter_behavior LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_NAME "test_muscl_reconstruction_limiter_behavior")

# Automatically find all source files in the current directory
AUX_SOURCE_DIRECTORY(. DIR_SRCS)

# Add the executable
ADD_EXECUTABLE(${PROJECT_NAME} ${DIR_SRCS})

# Link against the main SPHinXsys library (2D version)
target_link_libraries(${PROJECT_NAME} PRIVATE sphinxsys_2d)

# Find GTest
find_package(GTest REQUIRED)

# Link GTest targets
target_link_libraries(${PROJECT_NAME} PRIVATE GTest::gtest_main GTest::gtest)

# Add test to CTest
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
