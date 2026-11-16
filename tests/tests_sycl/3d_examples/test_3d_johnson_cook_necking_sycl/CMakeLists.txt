add_executable(test_3d_johnson_cook_necking_sycl necking_sycl.cpp)
target_link_libraries(test_3d_johnson_cook_necking_sycl PRIVATE sphinxsys_3d)
set(JC_CLASSIC_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/../../../3d_examples/test_3d_johnson_cook_necking")
target_include_directories(test_3d_johnson_cook_necking_sycl PRIVATE "${JC_CLASSIC_EXAMPLE}")
set_target_properties(test_3d_johnson_cook_necking_sycl PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/smoke")
if(SPHINXSYS_USE_FLOAT)
    set(JC_REGRESSION_PRECISION float)
else()
    set(JC_REGRESSION_PRECISION double)
endif()
file(COPY "${JC_CLASSIC_EXAMPLE}/regression_test_tool/${JC_REGRESSION_PRECISION}/"
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/smoke/input")
add_test(NAME test_3d_johnson_cook_necking_sycl
    COMMAND test_3d_johnson_cook_necking_sycl --regression-test --spacing=0.0002 --end-time=0.00008 --frames=4 --threads=2
        --state_recording=$<IF:$<BOOL:${TEST_STATE_RECORDING}>,true,false>)
set_tests_properties(test_3d_johnson_cook_necking_sycl PROPERTIES
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/smoke"
    LABELS "solid_dynamics;plasticity;computing_kernel" TIMEOUT 300)
