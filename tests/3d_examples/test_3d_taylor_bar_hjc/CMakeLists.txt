add_executable(test_3d_taylor_bar_hjc taylor_bar_hjc.cpp)
target_link_libraries(test_3d_taylor_bar_hjc PRIVATE sphinxsys_3d)
set_target_properties(test_3d_taylor_bar_hjc PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/smoke")
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/regression_test_tool/ DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/smoke/input")
add_test(NAME test_3d_taylor_bar_hjc COMMAND test_3d_taylor_bar_hjc --regression-test --spacing=0.001 --end-time=0.00006 --state_recording=false)
set_tests_properties(test_3d_taylor_bar_hjc PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/smoke"
    LABELS "solid_dynamics;contact" TIMEOUT 300)
