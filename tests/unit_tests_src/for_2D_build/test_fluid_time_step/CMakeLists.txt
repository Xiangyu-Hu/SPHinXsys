cmake_minimum_required(VERSION 3.10)

if(SPHINXSYS_2D AND SPHINXSYS_BUILD_2D_EXAMPLES)
    find_package(GTest REQUIRED)
    add_executable(test_fluid_time_step test_fluid_time_step.cpp)
    target_link_libraries(test_fluid_time_step PRIVATE sphinxsys_2d GTest::gtest_main GTest::gtest)
    add_test(NAME FluidTimeStepCompatibility COMMAND test_fluid_time_step)
endif()
