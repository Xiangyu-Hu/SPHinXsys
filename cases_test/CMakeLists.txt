SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})

#FOREACH(subdir ${SUBDIRS})
#	#MESSAGE("${subdir}")
#	ADD_SUBDIRECTORY(${subdir})
#ENDFOREACH()

ADD_SUBDIRECTORY(test_3d_muscle_compression)
ADD_SUBDIRECTORY(test_3d_particle_generation)