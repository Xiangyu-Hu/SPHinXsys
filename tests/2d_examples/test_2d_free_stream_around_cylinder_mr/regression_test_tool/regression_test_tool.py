# !/usr/bin/env python3
import os
import sys

path = os.path.abspath('../../../../../PythonScriptStore/RegressionTest')
sys.path.append(path)
from regression_test_base_tool import SphinxsysRegressionTest

"""
case name: test_2d_free_stream_around_cylinder_mr
"""

case_name = "test_2d_free_stream_around_cylinder_mr"
body_name = "Cylinder"
parameter_name = "TotalViscousForceFromFluid"

number_of_run_times = 0
converged = 0
sphinxsys = SphinxsysRegressionTest(case_name, body_name, parameter_name)

while True:
    print("Now start a new run......")
    sphinxsys.run_particle_relaxation()
    sphinxsys.run_case_with_reload()
    number_of_run_times += 1
    converged = sphinxsys.read_dat_file()
    print("Please note: This is the run number", number_of_run_times,".")
    if number_of_run_times <= 200:
        if converged == "true":
            print("The tested parameters of all variables are converged, and the run will stop here!")
            break
        elif converged != "true":
            print("The tested parameters of", sphinxsys.sphinxsys_parameter_name, "are not converged!")
            continue
    else:
        print("It's too many runs but still not converged, please try again!")
        break
