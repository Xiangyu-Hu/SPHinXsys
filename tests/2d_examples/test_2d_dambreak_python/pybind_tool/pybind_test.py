#!/usr/bin/env python3
import os
import sys
import platform
import argparse
# add dynamic link library or shared object to python env
# attention: match current python version with the version exposing the cpp code
sys_str = platform.system()
path_1 = os.path.abspath(os.path.join(os.getcwd(), '..'))
if sys_str == 'Windows':
    path_2 = 'lib'
elif sys_str == 'Linux':
    path_2 = 'lib'
else:
    # depend on the system
    path_2 = 'lib'
path = os.path.join(path_1, path_2)
sys.path.append(path)
# change import depending on the project name
import test_2d_dambreak_python as test_2d


def run_case():
    parser = argparse.ArgumentParser()
    # set case parameters
    parser.add_argument("--restart_step", default=0, type=int)
    parser.add_argument("--end_time", default=20, type=int)
    case = parser.parse_args()
    
    # set project from class, which is set in cpp pybind module
    project = test_2d.dambreak_from_sph_cpp(case.restart_step)
    if project.CmakeTest() == 1:
        project.RunCase(case.end_time)
    else:
        print("check path: ", path)
        

if __name__ == "__main__":
    run_case()
