STRING(REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})
PROJECT("${CURRENT_FOLDER}")

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")

add_executable(${PROJECT_NAME} test_ophelie_wcsph_operator_validation.cpp)
# Do NOT link extra_sources_2d: it currently pulls 3D-only french geometry into 2D.
target_link_libraries(${PROJECT_NAME} sphinxsys_2d)
target_include_directories(${PROJECT_NAME} PRIVATE
    ${SPHINXSYS_PROJECT_DIR}/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/hydrodynamics
    ${SPHINXSYS_PROJECT_DIR}/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie)

if(SPHINXSYS_USE_SYCL)
    file(GLOB UR_LOADER_CANDIDATES "/opt/intel/oneapi/compiler/*/lib/libur_loader.so*")
    list(SORT UR_LOADER_CANDIDATES)
    if(UR_LOADER_CANDIDATES)
        set(UR_LOADER_LIBRARY "")
        foreach(_cand IN LISTS UR_LOADER_CANDIDATES)
            if(_cand MATCHES "libur_loader\\.so\\.0\\.11")
                set(UR_LOADER_LIBRARY "${_cand}")
                break()
            endif()
        endforeach()
        if(NOT UR_LOADER_LIBRARY)
            list(GET UR_LOADER_CANDIDATES -1 UR_LOADER_LIBRARY)
        endif()
        target_link_libraries(${PROJECT_NAME} ${UR_LOADER_LIBRARY})
    endif()
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")

add_test(NAME ${PROJECT_NAME}_o0a COMMAND ${PROJECT_NAME} --case o0a
    --out ${EXECUTABLE_OUTPUT_PATH}/o0a_host_report.json
    WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
add_test(NAME ${PROJECT_NAME}_o0b COMMAND ${PROJECT_NAME} --case o0b
    --out ${EXECUTABLE_OUTPUT_PATH}/o0b_ck_report.json
    WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
set_tests_properties(${PROJECT_NAME}_o0a ${PROJECT_NAME}_o0b PROPERTIES LABELS "ophelie;doc203;o0")
