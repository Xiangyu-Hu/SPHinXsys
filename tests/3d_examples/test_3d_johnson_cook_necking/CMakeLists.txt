add_executable(test_3d_johnson_cook_necking necking.cpp)
target_link_libraries(test_3d_johnson_cook_necking PRIVATE sphinxsys_3d)
set_target_properties(test_3d_johnson_cook_necking PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/smoke")
if(SPHINXSYS_USE_FLOAT)
    set(JC_REGRESSION_PRECISION float)
else()
    set(JC_REGRESSION_PRECISION double)
endif()
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/regression_test_tool/${JC_REGRESSION_PRECISION}/"
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/smoke/input")
add_test(NAME test_3d_johnson_cook_necking
    COMMAND test_3d_johnson_cook_necking --regression-test --spacing=0.0002 --end-time=0.00008 --frames=4 --threads=2
        --state_recording=$<IF:$<BOOL:${TEST_STATE_RECORDING}>,true,false>)
set_tests_properties(test_3d_johnson_cook_necking PROPERTIES
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/smoke"
    LABELS "solid_dynamics;plasticity" TIMEOUT 180)
