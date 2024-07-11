# install tianshou library
import os
import sys
import platform
import argparse


sys_str = platform.system()
# swith to current position install cd pip install -e .
#!/usr/bin/env python3

# add dynamic link library or shared object to python env
# attention: match current python version with the version exposing the cpp code

# If this doesn't works, try path_1 = os.path.abspath(os.path.join(os.getcwd(), '../..'))
def list_files_in_current_directory(path):
    # 获取当前目录
    # 列出当前目录中的所有文件和文件夹
    if os.path.exists(path):
        files_and_dirs = os.listdir(path)
        files = [f for f in files_and_dirs if os.path.isfile(os.path.join(path, f))]
        
        # 打印文件名
        for file in files:
            print(file)
        return True
    else:
        print(f"The directory {path} does not exist.")
        return False

path_1 = os.path.abspath(os.path.join(os.getcwd(), '../..'))
# path_1 = os.path.abspath(os.path.join(os.getcwd(), '..'))


print('path1:', path_1)
if sys_str == 'Windows':
    path_2 = 'Release'
elif sys_str == 'Linux':
    path_2 = ''
else:
    path_2 = 'lib'

path = os.path.join(path_1, path_2)
print(path)
sys.path.append(path)
# change import depending on the project name

list_files_in_current_directory(path_1)
if list_files_in_current_directory(path) == True:
    print("import test_2d_free_stream_around_fish_test_pybind")
    import test_2d_free_stream_around_fish_test_pybind as train
    r = train.from_sph_relaxation_and_test(2)
    b = train.from_sph_reload_and_test(2)
    b.SetFreq(4)
    b.RunCase(1, 0.1)
    print("success")
else:
    print("error")

# 在cmake文件中生成 一个是加测试 运行add_test  DQN文件
# 在cmake文件中
#add_test(NAME ${PROJECT_NAME} COMMAND  ${Python3_EXECUTABLE} "${EXECUTABLE_OUTPUT_PATH}/bind/pybind_test.py")
#set_tests_properties(${PROJECT_NAME} PROPERTIES WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}"
#    PASS_REGULAR_EXPRESSION "The result of Pressure is correct based on the dynamic time warping regression test!")