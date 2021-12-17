# !/usr/bin/env python3
import os


class SphinxsysRegressionTest:

    def __init__(self, casename, bodyname, parametername):
        self.sphinxsys_exec_path = os.path.abspath(os.path.join(os.getcwd()))
        self.sphinxsys_case_path = os.path.abspath(os.path.join(self.sphinxsys_exec_path))
        self.sphinxsys_case_name = casename
        self.sphinxsys_body_name = bodyname
        self.sphinxsys_parameter_name = parametername
        self.enter_sphinxsys_case_folder = f"cd {self.sphinxsys_case_path};"
        self.input_file_path = os.path.join(self.sphinxsys_case_path, "input")
        self.condition_file_path = os.path.join(self.input_file_path, f"{bodyname}_{parametername}_runtimes.dat")

    def compile_case(self) -> None:
        print('Start compiling test case....')
        command = "make -j8"
        os.system(self.enter_sphinxsys_case_folder + command)
        print('Compiling test case is finished...')

    def run_case(self) -> None:
        print('Start case simulation...')
        command = f"./{self.sphinxsys_case_name}"
        os.system(self.enter_sphinxsys_case_folder + command)
        print('Simulating case is finished...')

    def read_dat_file(self):
        file = open(self.condition_file_path)
        ifconverged = file.readline(4)
        file.close()
        return ifconverged


"""
case name: test_0d_regression_test
"""

case_name = "test_0d_regression_test"
body_name = "DiffusionBody"
parameter_name = "TotalAveragedParameterOnPartlyDiffusionBody"
body_name_1 = "TemperatureObserver"
parameter_name_1 = "Phi"

number_of_run_times = 0
converged = 0
sphinxsys = SphinxsysRegressionTest(case_name, body_name, parameter_name)
sphinxsys_1 = SphinxsysRegressionTest(case_name, body_name_1, parameter_name_1)
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
