# modules/Palace/CMakeLists.txt

# If SPHinXsys is not built with Palace support, exit early.
# This allows optional compilation of the Palace module.
if (NOT SPHINXSYS_WITH_PALACE)
    return()
endif()

# Define interface library for Palace integration.
# INTERFACE means this target only carries usage requirements
# (include dirs, link libs) but does not compile source itself. （Temporary）
add_library(sphinxsys_palace INTERFACE
    #palace_interface.cpp
)

# Link Palace::palace (defined by FindPalace.cmake) so any target
# linking to sphinxsys_palace will automatically get all Palace libs.  （Temporary）
target_link_libraries(sphinxsys_palace
    INTERFACE
    #PUBLIC
        Palace::palace
)
# Export include directory to users of sphinxsys_palace.
target_include_directories(sphinxsys_palace
    INTERFACE
    #PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Add test subdirectory that contains test cases for Palace.
add_subdirectory(tests)

