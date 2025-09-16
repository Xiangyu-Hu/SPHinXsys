# MUSCL reconstruction unit test
find_package(GTest REQUIRED)
add_executable(test_muscl_reconstruction
    test_muscl_reconstruction_linear.cpp
)

target_link_libraries(test_muscl_reconstruction
    PRIVATE sphinxsys_2d
    PRIVATE GTest::gtest_main
    PRIVATE GTest::gtest
)

add_test(NAME MUSCLReconstructionTest
    COMMAND test_muscl_reconstruction
)
