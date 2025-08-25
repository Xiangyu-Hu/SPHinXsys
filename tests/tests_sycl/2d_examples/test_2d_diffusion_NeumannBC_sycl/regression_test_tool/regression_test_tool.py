# !/usr/bin/env python3
import os
import sys

def find_upper_folder(folder_name, start_path=None):
    """
    Search for an upper-level folder by name, starting from a given directory and traversing upwards.
    """
    if start_path is None:
        start_path = os.getcwd()

    current_path = start_path
    while True:
        potential_path = os.path.join(current_path, folder_name)
        if os.path.isdir(potential_path):
            return potential_path
        parent_path = os.path.dirname(current_path)
        if parent_path == current_path:  # Reached the root directory
            break
        current_path = parent_path

    return None

def find_nested_folder(folder_chain, start_path=None):
    """
    Find a nested folder by traversing upward and then descending into a specified folder chain.
    """
    if not folder_chain:
        return None

    # Step 1: Find the first folder in the chain
    first_folder = folder_chain[0]
    first_folder_path = find_upper_folder(first_folder, start_path)
    if not first_folder_path:
        return None

    # Step 2: Traverse the rest of the folder chain
    current_path = first_folder_path
    for folder in folder_chain[1:]:
        current_path = os.path.join(current_path, folder)
        if not os.path.isdir(current_path):
            return None

    return current_path

# Step 1: Find the RegressionTest folder
folder_chain = ["build", "PythonScriptStore", "RegressionTest"]
regression_test_folder = find_nested_folder(folder_chain)

if not regression_test_folder:
    print("'RegressionTest' folder not found.")
    exit(1)

# Step 2: Add the RegressionTest folder to sys.path
if regression_test_folder not in sys.path:
    sys.path.insert(0, regression_test_folder)

# Step 3: Import the module
try:
    from regression_test_base_tool import SphinxsysRegressionTest
    print("Module imported successfully!")
except ImportError as e:
    print(f"Failed to import module: {e}")


"""
case name: test_2d_diffusion_NeumannBC
"""

case_name = "test_2d_diffusion_NeumannBC_sycl"
body_name = "TemperatureObserver"
parameter_name = "Phi"

number_of_run_times = 0
converged = 0
sphinxsys = SphinxsysRegressionTest(case_name, body_name, parameter_name)


while True:
    print("Now start a new run......")
    sphinxsys.run_case()
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
        print("It's too many runs but still not converged, please try again!")
        break
