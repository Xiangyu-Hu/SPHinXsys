add_executable(test_3d_hjc_impact hjc_impact.cpp)
target_link_libraries(test_3d_hjc_impact PRIVATE sphinxsys_3d)
set_target_properties(test_3d_hjc_impact PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/smoke")
add_test(NAME test_3d_hjc_impact COMMAND test_3d_hjc_impact --spacing=0.001 --end-time=0.00006 --state_recording=false)
set_tests_properties(test_3d_hjc_impact PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/smoke"
    LABELS "solid_dynamics;contact" TIMEOUT 300)
