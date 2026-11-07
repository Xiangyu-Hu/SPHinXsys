STRING(REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})
PROJECT("${CURRENT_FOLDER}")

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")

add_executable(${PROJECT_NAME} test_2d_ophelie_natural_convection_validation.cpp)
target_link_libraries(${PROJECT_NAME} sphinxsys_2d)
target_include_directories(${PROJECT_NAME} PRIVATE
    ${SPHINXSYS_PROJECT_DIR}/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie
    ${SPHINXSYS_PROJECT_DIR}/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/diagnostics
    ${SPHINXSYS_PROJECT_DIR}/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/hydrodynamics)

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
