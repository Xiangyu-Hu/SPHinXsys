# install tianshou library
import os
import sys
import subprocess
import platform
import argparse


sys_str = platform.system()
# def install(package):
#     subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# def install_self():
#     subprocess.check_call([sys.executable, "-m", "pip", "install", "-e", "."])
# try:
#     import tianshou
# except ImportError:
#     if sys_str == 'Windows':
#         print("tianshou library not found. Installing...")
#         install("tianshou")
#         import tianshou
#     elif sys_str == 'Linux':
#         print("numpy library reinstalling...")
#         install("numpy==1.23.4")
#         print("tianshou library not found. Installing...")
#         install("tianshou")
#         print("finish install numpy & tianshou")
#         try:
#             import tianshou
#         except ImportError:
#             print("install tianshou failed ...")
#     else:
#         print("tianshou library not found. Installing...")
#         install("tianshou")

# swith to current position install cd pip install -e .
#!/usr/bin/env python3

# add dynamic link library or shared object to python env
# attention: match current python version with the version exposing the cpp code

# If this doesn't works, try path_1 = os.path.abspath(os.path.join(os.getcwd(), '../..'))
def list_files_in_current_directory(path):
    # 获取当前目录
    current_directory = os.getcwd()
    # 列出当前目录中的所有文件和文件夹
    files_and_dirs = os.listdir(path)
    
    # 过滤出文件
    files = [f for f in files_and_dirs if os.path.isfile(os.path.join(path, f))]
    
    # 打印文件名
    for file in files:
        print(file)


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
# change import depending on the project name
list_files_in_current_directory(path)

print("import test_2d_free_stream_around_fish_pybind")
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