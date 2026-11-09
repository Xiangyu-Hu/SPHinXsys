SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})

foreach(subdir ${SUBDIRS})
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt)
        # With SYCL, only GPU-paper instrumented cases (legacy extras may be non-SYCL).
        if(SPHINXSYS_USE_SYCL AND NOT subdir MATCHES "gpu_paper")
            continue()
        endif()
        add_subdirectory(${subdir})
    endif()
endforeach()
