cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

if(SPHINXSYS_2D AND SPHINXSYS_BUILD_2D_EXAMPLES)

project(test_muscl_reconstruction_limiter_behavior LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_NAME "test_muscl_reconstruction_limiter_behavior")

AUX_SOURCE_DIRECTORY(. DIR_SRCS)

add_executable(${PROJECT_NAME} ${DIR_SRCS})

find_package(GTest REQUIRED)

target_link_libraries(${PROJECT_NAME}
    PRIVATE sphinxsys_2d
    PRIVATE GTest::gtest_main
    PRIVATE GTest::gtest
)

add_test(
    NAME ${PROJECT_NAME}
    COMMAND ${PROJECT_NAME}
)

endif()