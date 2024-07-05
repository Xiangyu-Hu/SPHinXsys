# install tianshou library
import os
import sys
import subprocess

def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

def install_self():
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-e", "."])
try:
    import tianshou
except ImportError:
    print("tianshou library not found. Installing...")
    install("tianshou")
    import tianshou

print("tianshou library is installed and ready to use.")

# swith to current position install cd pip install -e .
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
print('path1:', path_1)
if sys_str == 'Windows':
    # Append 'RelWithDebInfo' or 'Debug' depending on the configuration
    # For example, path_2 = 'lib/RelWithDebInfo'
    path_2 = 'Release'
elif sys_str == 'Linux':
    path_2 = 'lib'
else:
    # depend on the system
    path_2 = 'lib'

path = os.path.join(path_1, path_2)
print(path)
sys.path.append(path)
# change import depending on the project name
# 切换到目标文件夹
path_train = 'gym-flap_fish'
path_test = 'gym-flap_fish_test'
target_train_directory = os.path.join(path_1, path_train)
target_test_directory = os.path.join(path_1, path_test)


print("gym_fish library Installing...")
target_train_directory = os.path.join(path_1, path_train)
print('target_train_directory: ', target_train_directory)
os.chdir(target_train_directory)
install_self()

print("flap_fish_test library Installing...")
target_test_directory = os.path.join(path_1, path_test)
print('target_test_directory: ', target_test_directory)
os.chdir(target_test_directory)
install_self()

# pip install -e 
# import test_2d_flow_stream_around_fish_flap_pybind as train
print("swith to run test algorithms!")
path_test = 'drl_tool'
target_directory = os.path.join(path_1, path_test)
print('target_directory: ', target_directory)
os.chdir(target_directory)
subprocess.check_call([sys.executable, "sac.py"])
# r = train.from_sph_relaxation_and_test(2)
# b = train.from_sph_reload_and_test(2)
# b.SetAmUp(0.12)
# b.SetAmDown(0.12)
# b.RunCase(1, 0.2)

# 在cmake文件中生成 一个是加测试 运行add_test  DQN文件
# 在cmake文件中
#add_test(NAME ${PROJECT_NAME} COMMAND  ${Python3_EXECUTABLE} "${EXECUTABLE_OUTPUT_PATH}/bind/pybind_test.py")
#set_tests_properties(${PROJECT_NAME} PROPERTIES WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}"
#    PASS_REGULAR_EXPRESSION "The result of Pressure is correct based on the dynamic time warping regression test!")