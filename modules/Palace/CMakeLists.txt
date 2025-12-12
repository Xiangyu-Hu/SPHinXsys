if (NOT SPHINXSYS_WITH_PALACE)
    return()
endif()

add_library(sphinxsys_palace
    sphinxsys_palace_interface.cpp
    palace_magnetostatic_interface.cpp
)

target_link_libraries(sphinxsys_palace
    PUBLIC
        Palace::palace
)

target_include_directories(sphinxsys_palace
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

add_subdirectory(tests)
