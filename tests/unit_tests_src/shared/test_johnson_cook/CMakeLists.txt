add_executable(test_johnson_cook test_johnson_cook.cpp)
target_link_libraries(test_johnson_cook PRIVATE sphinxsys_3d GTest::gtest_main)
gtest_discover_tests(test_johnson_cook TEST_PREFIX "JohnsonCook.")
