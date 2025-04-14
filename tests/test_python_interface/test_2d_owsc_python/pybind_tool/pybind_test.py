#!/usr/bin/env python3
import os
import sys
import platform
import argparse
# add dynamic link library or shared object to python env
# attention: match current python version with the version exposing the cpp code
sys_str = platform.system()
# If this doesn't works, try path_1 = os.path.abspath(os.path.join(os.getcwd(), '../..'))
path_1 = os.path.abspath(os.path.join(os.getcwd(), '..'))
if sys_str == 'Windows':
    # Append 'RelWithDebInfo' or 'Debug' depending on the configuration
    # For example, path_2 = 'lib/RelWithDebInfo'
    path_2 = 'lib'
elif sys_str == 'Linux':
    path_2 = 'lib'
else:
    # depend on the system
    path_2 = 'lib'
path = os.path.join(path_1, path_2)
sys.path.append(path)
# change import depending on the project name
import test_2d_owsc_python as test_2d


def run_case():
    parser = argparse.ArgumentParser()
    # set case parameters
    parser.add_argument("--parallel_env", default=0, type=int)
    parser.add_argument("--episode_env", default=0, type=int)
    parser.add_argument("--damping_coefficient", default=20.0, type=float)
    parser.add_argument("--end_time", default=12.0, type=float)
    case = parser.parse_args()
    
    # set project from class, which is set in cpp pybind module
    project = test_2d.owsc_from_sph_cpp(case.parallel_env, case.episode_env)
    if project.cmake_test() == 1:
        project.run_case(case.end_time, case.damping_coefficient)
    else:
        print("check path: ", path)
        

if __name__ == "__main__":
    run_case()
