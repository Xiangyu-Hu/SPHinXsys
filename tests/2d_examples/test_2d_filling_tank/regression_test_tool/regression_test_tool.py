# !/usr/bin/env python3
import os


class SphinxsysRegressionTest:

    def __init__(self, casename, bodyname, parametername):
        self.sphinxsys_exec_path = os.path.abspath(os.path.join(os.getcwd(), ".."))
        self.sphinxsys_case_path = os.path.abspath(os.path.join(self.sphinxsys_exec_path, ".."))
        self.sphinxsys_src_path = os.path.join(self.sphinxsys_case_path, "src")
        self.sphinxsys_rld_path = os.path.join(self.sphinxsys_src_path, "reload")
        self.sphinxsys_case_name = casename
        self.sphinxsys_body_name = bodyname
        self.sphinxsys_parameter_name = parametername
        self.enter_sphinxsys_exec_folder = f"cd {self.sphinxsys_exec_path};"
        self.enter_sphinxsys_case_folder = f"cd {self.sphinxsys_case_path};"
        self.input_file_path = os.path.join(self.sphinxsys_exec_path, "input")
        self.condition_file_path = os.path.join(self.input_file_path, f"{bodyname}_{parametername}_runtimes.dat")

    def compile_case(self) -> None:
        print('Start compiling test case....')
        command = "make -j8"
        os.system(self.enter_sphinxsys_case_folder + command)
        print('Compiling test case is finished...')

    def test_case(self) -> None:
        print('Start test case...')
        command = "make test"
        os.system(self.enter_sphinxsys_case_folder + command)
        print('Testing case is finished...')

    def copy_reload(self) -> None:
        print('Start copy the reload file...')
        command = "cp -r reload bin"
        os.system(self.enter_sphinxsys_case_folder + command)
        print('Copying threload file is finished...')

    def run_case(self) -> None:
        print('Start case simulation...')
        command = f"./{self.sphinxsys_case_name}"
        os.system(self.enter_sphinxsys_exec_folder + command)
        print('Simulating case is finished...')

    def read_dat_file(self):
        file = open(self.condition_file_path)
        ifconverged = file.readline(4)
        file.close()
        return ifconverged


"""
case name: test_2d_filling_tank
"""

case_name = "test_2d_filling_tank"
body_name = "WaterBody"
parameter_name = "TotalMechanicalEnergy"
body_name_1 = "Fluidobserver"
parameter_name_1 = "Pressure"

number_of_run_times = 0
converged = 0
sphinxsys = SphinxsysRegressionTest(case_name, body_name, parameter_name)
sphinxsys_1 = SphinxsysRegressionTest(case_name, body_name_1, parameter_name_1)


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
