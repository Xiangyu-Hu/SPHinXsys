project(test_3d_taylor_bar_hjc_sycl)
add_executable(${PROJECT_NAME} taylor_bar_hjc_sycl.cpp)
target_link_libraries(${PROJECT_NAME} sphinxsys_3d)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/input)
# Reuse the classic example's baseline for both precisions. Only the mean-damage
# recording name differs; its numerical means, variances and run counts stay unchanged.
set(HJC_CLASSIC_REFERENCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../3d_examples/test_3d_taylor_bar_hjc/regression_test_tool")
foreach(HJC_REFERENCE_SUFFIX ensemble_averaged_mean_variance.xml runtimes.dat)
    configure_file("${HJC_CLASSIC_REFERENCE_DIR}/Concrete_TotalKineticEnergy_${HJC_REFERENCE_SUFFIX}"
                   "${CMAKE_CURRENT_BINARY_DIR}/input/Concrete_TotalKineticEnergy_${HJC_REFERENCE_SUFFIX}" COPYONLY)
    set(HJC_DAMAGE_INPUT "${CMAKE_CURRENT_BINARY_DIR}/input/Concrete_AverageHJCDamage_${HJC_REFERENCE_SUFFIX}")
    configure_file("${HJC_CLASSIC_REFERENCE_DIR}/Concrete_TotalHJCDamage_${HJC_REFERENCE_SUFFIX}"
                   "${HJC_DAMAGE_INPUT}" COPYONLY)
    file(READ "${HJC_DAMAGE_INPUT}" HJC_DAMAGE_REFERENCE)
    string(REPLACE "TotalHJCDamage" "AverageHJCDamage" HJC_DAMAGE_REFERENCE "${HJC_DAMAGE_REFERENCE}")
    file(WRITE "${HJC_DAMAGE_INPUT}" "${HJC_DAMAGE_REFERENCE}")
endforeach()
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} --spacing=0.001 --regression-test --state_recording=${TEST_STATE_RECORDING}
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(${PROJECT_NAME} PROPERTIES LABELS "solid_dynamics; contact")
