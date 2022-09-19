# !/usr/bin/env python3
import os
import sys

path = os.path.abspath('../../../../PythonScriptStore/RegressionTest')
sys.path.append(path)
print(path)
from regression_test_base_tool import SphinxsysRegressionTestByCTest

"""
case name: test_0d_regression_test
"""

case_name = "test_0d_regression_test"
body_name = "InnerDomain"
parameter_name = "DiffusionReactionSpeciesAverage"
body_name_1 = "TemperatureObserver"
parameter_name_1 = "Phi"

number_of_run_times = 0
converged = 0
sphinxsys = SphinxsysRegressionTestByCTest(case_name, body_name, parameter_name)
sphinxsys_1 = SphinxsysRegressionTestByCTest(case_name, body_name_1, parameter_name_1)
sphinxsys.run_case()


while True:
    print("Now start a new run......")
    sphinxsys.run_case()
    number_of_run_times += 1
    converged = sphinxsys.read_dat_file()
    converged_1 = sphinxsys_1.read_dat_file()
    print("Please note: This is the", number_of_run_times, "run!")
    if number_of_run_times <= 200:
        if (converged == "true") and (converged_1 == "true"):
            print("The tested parameters of all variables are converged, and the run will stop here!")
            break
        elif converged != "true":
            print("The tested parameters of", sphinxsys.sphinxsys_parameter_name, "are not converged!")
            continue
        elif converged_1 != "true":
            print("The tested parameters of", sphinxsys_1.sphinxsys_parameter_name, "are not converged!")
            continue
    else:
        print("It's too many runs but still not converged, please try again!")
        break
