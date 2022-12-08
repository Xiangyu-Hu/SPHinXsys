
if(SPHINXSYS_2D AND SPHINXSYS_BUILD_2D_EXAMPLES)
	ADD_SUBDIRECTORY(for_2D_build)
endif()

if(SPHINXSYS_3D AND SPHINXSYS_BUILD_3D_EXAMPLES)
	ADD_SUBDIRECTORY(shared) # TODO make a base target containing the shared code btw 2d and 3d 
	ADD_SUBDIRECTORY(for_3D_build)
endif()
