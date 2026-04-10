cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

if(SPHINXSYS_2D AND SPHINXSYS_BUILD_2D_EXAMPLES)

project(test_muscl_hllc_integration LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(test_muscl_hllc_integration
    ${CMAKE_CURRENT_SOURCE_DIR}/test_muscl_hllc_integration.cpp
)

find_package(GTest REQUIRED)

target_link_libraries(test_muscl_hllc_integration
    PRIVATE sphinxsys_2d
    PRIVATE GTest::gtest_main
    PRIVATE GTest::gtest
)

add_test(
    NAME test_muscl_hllc_integration
    COMMAND test_muscl_hllc_integration
)

endif()