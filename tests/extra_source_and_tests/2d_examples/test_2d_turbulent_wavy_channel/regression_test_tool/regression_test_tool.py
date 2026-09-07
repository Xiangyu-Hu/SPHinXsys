# !/usr/bin/env python3
import os
import sys

path = os.path.abspath('../../../../../../PythonScriptStore/RegressionTest')
sys.path.append(path)
from regression_test_base_tool import SphinxsysRegressionTest

"""
case name: test_2d_turbulent_wavy_channel
"""

case_name = "test_2d_turbulent_wavy_channel"
body_name = "ObserverCenterPoint"
parameter_name = "TurbulentViscosity"


number_of_run_times = 0
converged = 0
sphinxsys = SphinxsysRegressionTest(case_name, body_name, parameter_name)


if not os.path.isfile(os.path.join(sphinxsys.sphinxsys_exec_path, "reload", "Reload.xml")):
    sphinxsys.run_particle_relaxation()


while number_of_run_times < 200:
    print("Now start a new run......")
    sphinxsys.run_case_with_reload()
    number_of_run_times += 1
    converged = sphinxsys.read_dat_file()
    print("Please note: This is the", number_of_run_times, "run!")
    if number_of_run_times <= 200:
        if (converged == "true"):
            print("The tested parameters of all variables are converged, and the run will stop here!")
            break
        elif converged != "true":
            print("The tested parameters of", sphinxsys.sphinxsys_parameter_name, "are not converged!")
            continue
else:
    raise SystemExit("Database not converged after 200 runs; review the results and threshold.")
