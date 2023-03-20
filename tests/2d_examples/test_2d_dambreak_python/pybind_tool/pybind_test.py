# !/usr/bin/env python3
import os
import sys
import platform
import argparse
# add dynamic link library or shared object to python env
# attention: match current python version with the version exposing the cpp code
sys_str = platform.system()
path_1 = os.path.abspath(os.path.join(os.getcwd(), '../..'))
if sys_str == 'Windows':
    path_2 = 'lib\\RelWithDebInfo'
elif sys_str == 'Linux':
    path_2 = 'lib'
else:
    # depend on the system
    path_2 = 'lib'
path = os.path.join(path_1, path_2)
sys.path.append(path)
# change import depending on the project name
import test_2d_dambreak_python as test_2d

ctest = test_2d.dambreak_from_sph_cpp(0)
result = ctest.CmakeTest()

return result



