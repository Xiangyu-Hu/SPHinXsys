STRING(REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})
PROJECT("${CURRENT_FOLDER}")

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")

add_executable(${PROJECT_NAME} test_3d_ophelie.cpp)
target_link_libraries(${PROJECT_NAME} sphinxsys_3d extra_sources_3d)

# SYCL linking may require explicit UR loader linkage.
# Your build log shows missing `ur* @LIBUR_LOADER_0.11` symbols from `libsycl.so`,
# which is provided by `libur_loader.so`.
if(SPHINXSYS_USE_SYCL)
    file(GLOB UR_LOADER_CANDIDATES "/opt/intel/oneapi/compiler/*/lib/libur_loader.so*")
    list(SORT UR_LOADER_CANDIDATES)
    if(UR_LOADER_CANDIDATES)
        # `libsycl.so` 里缺失的符号带有 `@LIBUR_LOADER_0.11` 版本标记，
        # 因此这里优先选 0.11.x 的具体库文件，避免链接到 newer ABI (e.g. 0.12.x)。
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
        message(STATUS "Linking with ${UR_LOADER_LIBRARY}")
        target_link_libraries(${PROJECT_NAME} ${UR_LOADER_LIBRARY})
    else()
        message(WARNING "libur_loader.so not found under /opt/intel/oneapi/compiler/*/lib; SYCL link may fail.")
    endif()
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")
