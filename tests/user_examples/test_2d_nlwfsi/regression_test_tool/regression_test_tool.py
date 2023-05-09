# !/usr/bin/env python3
import os
import sys

path = os.path.abspath('../../../../../PythonScriptStore/RegressionTest')
sys.path.append(path)
from regression_test_base_tool import SphinxsysRegressionTest

"""
case name: test_2d_nlwfsi
"""

case_name = "test_2d_nlwfsi"
body_name = "WaveGauge"
parameter_name = "FreeSurfaceHeight"
body_name_1 = "Observer"
parameter_name_1 = "Position"
body_name_2 = "FluidObserver2"
parameter_name_2 = "Pressure"
body_name_3 = "FluidObserver3"
parameter_name_3 = "Pressure"

number_of_run_times = 0
converged = 0
sphinxsys = SphinxsysRegressionTest(case_name, body_name, parameter_name)
sphinxsys_1 = SphinxsysRegressionTest(case_name, body_name_1, parameter_name_1)
sphinxsys_2 = SphinxsysRegressionTest(case_name, body_name_2, parameter_name_2)
sphinxsys_3 = SphinxsysRegressionTest(case_name, body_name_3, parameter_name_3)

while True:
    print("Now start a new run......")
    sphinxsys.run_case()
    number_of_run_times += 1
    converged = sphinxsys.read_dat_file()
    converged_1 = sphinxsys_1.read_dat_file()
    converged_2 = sphinxsys_2.read_dat_file()
    converged_3 = sphinxsys_3.read_dat_file()
    print("Please note: This is the", number_of_run_times, "run!")
    if number_of_run_times <= 200:
        if (converged == "true") and (converged_1 == "true") and (converged_2 == "true") and (converged_3 == "true"):
            print("The tested parameters of all variables are converged, and the run will stop here!")
            break
        elif converged != "true":
            print("The tested parameters of", sphinxsys.sphinxsys_parameter_name, "are not converged!")
            continue
        elif converged_1 != "true":
            print("The tested parameters of", sphinxsys_1.sphinxsys_parameter_name, "are not converged!")
            continue
        elif converged_2 != "true":
            print("The tested parameters of", sphinxsys_2.sphinxsys_parameter_name, "are not converged!")
            continue
        elif converged_3 != "true":
            print("The tested parameters of", sphinxsys_3.sphinxsys_parameter_name, "are not converged!")
            continue
    else:
        print("It's too many runs but still not converged, please try again!")
        break
