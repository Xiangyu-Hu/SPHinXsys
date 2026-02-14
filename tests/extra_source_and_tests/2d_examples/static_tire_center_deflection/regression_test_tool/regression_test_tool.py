
import os
import sys

path = os.path.abspath('../../../../../PythonScriptStore/RegressionTest')
sys.path.append(path)

from regression_test_base_tool import SphinxsysRegressionTest


case_name = "static_tire_center_deflection"
body_name = "TireObserver"
parameter_name = "Position"
converged = 0
number_of_run_times = 0
sphinxsys = SphinxsysRegressionTest(case_name, body_name, parameter_name)


clean_input_folder(sphinxsys.input_file_path)

while True:
    print("Now start a new run......")
    sphinxsys.run_case()
    number_of_run_times += 1
    converged = sphinxsys.read_dat_file()
    print("Please note: This is the", number_of_run_times, "run!")

    if number_of_run_times <= 200:
        if converged == "true":
            print("Converged. Stop.")
            break
        elif converged != "true":
            print("Not converged, continue...")
            continue
    else:
        print("Too many runs, still not converged.")
        break


