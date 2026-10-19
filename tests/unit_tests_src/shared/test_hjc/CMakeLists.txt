add_executable(test_hjc test_hjc.cpp)
target_link_libraries(test_hjc PRIVATE sphinxsys_3d GTest::gtest_main)
gtest_discover_tests(test_hjc TEST_PREFIX "HJC.")
