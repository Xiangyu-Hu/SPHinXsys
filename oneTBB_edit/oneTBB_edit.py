# this is a python file to adjust oneTBB
import shutil

##### oneTBB main CmakeLists changes #####

# turn off testing: option(TBB_TEST "Enable testing" OFF)
shutil.copy('CMakeLists.txt', '../oneTBB')

##### oneTBB main CmakeLists changes #####

##### tbbmalloc.cpp changes #####

# remove: #include "../tbb/assert_impl.h"
shutil.copy('tbbmalloc.cpp', '../oneTBB/src/tbbmalloc')

##### tbbmalloc.cpp changes #####