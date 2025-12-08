# modules/Palace/CMakeLists.txt


if (NOT SPHINXSYS_WITH_PALACE)
    return()
endif()


add_library(sphinxsys_palace INTERFACE
    #palace_interface.cpp
)


target_link_libraries(sphinxsys_palace
    INTERFACE
    #PUBLIC
        Palace::palace
)

target_include_directories(sphinxsys_palace
    INTERFACE
    #PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)


add_subdirectory(tests)

