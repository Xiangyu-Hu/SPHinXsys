if(SPHINXSYS_2D)
    add_executable(test_shepard_density_regularization
        test_shepard_density_regularization.cpp)
    target_link_libraries(test_shepard_density_regularization
        PRIVATE sphinxsys_2d GTest::gtest_main GTest::gtest)
    add_test(NAME test_shepard_density_regularization
        COMMAND test_shepard_density_regularization)
endif()
