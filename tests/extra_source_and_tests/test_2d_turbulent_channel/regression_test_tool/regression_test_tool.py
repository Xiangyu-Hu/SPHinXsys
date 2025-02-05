# !/usr/bin/env python3
import os
import sys

path = os.path.abspath('../../../../../PythonScriptStore/RegressionTest')
sys.path.append(path)
from regression_test_base_tool import SphinxsysRegressionTest

"""
case name: 2d_turbulent_channel_VI_offset_model
"""

case_name = "2d_turbulent_channel_VI_offset_model"
body_name = "ObserverCenterPoint"
parameter_name = "TurbulentViscosity"


number_of_run_times = 0
converged = 0
sphinxsys = SphinxsysRegressionTest(case_name, body_name, parameter_name)


while True:
    print("Now start a new run......")
    sphinxsys.run_case()
    number_of_run_times += 1
    converged = sphinxsys.read_dat_file()
    print("Please note: This is the", number_of_run_times, "run!")
    if number_of_run_times <= 100:
        if (converged == "true"):
            print("The tested parameters of all variables are converged, and the run will stop here!")
            break
        elif converged != "true":
            print("The tested parameters of", sphinxsys.sphinxsys_parameter_name, "are not converged!")
            continue
    else:
        print("It's too many runs but still not converged, please try again!")
        break
