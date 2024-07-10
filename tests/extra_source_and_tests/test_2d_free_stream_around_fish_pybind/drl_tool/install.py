# install tianshou library
import os
import sys
import subprocess
import platform
import argparse


sys_str = platform.system()

# swith to current position install cd pip install -e .
#!/usr/bin/env python3

# add dynamic link library or shared object to python env
# attention: match current python version with the version exposing the cpp code

# If this doesn't works, try path_1 = os.path.abspath(os.path.join(os.getcwd(), '../..'))
def read_files_in_folder(folder_path):
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            file_path = os.path.join(root, file)
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    print(f"Contents of {file_path}:")
                    print(f.read())
                    print("=" * 80)  # Separator between file contents
            except Exception as e:
                print(f"Could not read {file_path} due to {e}")
path_1 = os.path.abspath(os.path.join(os.getcwd(), '..'))


print('path1:', path_1)
if sys_str == 'Windows':
    path_2 = 'Release'
elif sys_str == 'Linux':
    path_2 = 'lib'
else:
    path_2 = 'lib'

path = os.path.join(path_1, path_2)
print(path)
sys.path.append(path)
read_files_in_folder(path)
# change import depending on the project name

import test_2d_free_stream_around_fish_pybind as train
r = train.from_sph_relaxation(2)
b = train.from_sph_reload_and_train(2)
b.SetFreq(4)
b.RunCase(1, 0.1)
print("success")


# 在cmake文件中生成 一个是加测试 运行add_test  DQN文件
# 在cmake文件中
#add_test(NAME ${PROJECT_NAME} COMMAND  ${Python3_EXECUTABLE} "${EXECUTABLE_OUTPUT_PATH}/bind/pybind_test.py")
#set_tests_properties(${PROJECT_NAME} PROPERTIES WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}"
#    PASS_REGULAR_EXPRESSION "The result of Pressure is correct based on the dynamic time warping regression test!")