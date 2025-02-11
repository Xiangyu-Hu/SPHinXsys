# !/usr/bin/env python3
import os

# This is the regression test base module used for all test cases.
# There are two different schemes that differ from the executable path.
# The first one is calling Python script by CTest, such as test_0d_regression_test.
# The other one is calling Python script directly manually.
# They have different relative paths, and correspond to different modules.


class SphinxsysRegressionTestByCTest:

    def __init__(self, casename, bodyname, parametername):
        self.sphinxsys_exec_path = os.path.abspath(os.path.join(os.getcwd()))
        self.sphinxsys_case_path = os.path.abspath(os.path.join(self.sphinxsys_exec_path))
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
        os.system(self.enter_sphinxsys_case_folder)
        os.system(command)
        print('Compiling test case is finished...')

    def run_particle_relaxation(self) -> None:
        print('Start particle relaxation for the simulation...')
        command = f".{os.sep}{self.sphinxsys_case_name} --relax=true"
        os.system(self.enter_sphinxsys_exec_folder)
        os.system(command)
        print('Simulating case is finished...')

    def run_case(self) -> None:
        print('Start case simulation...')
        print(self.enter_sphinxsys_exec_folder)
        command = f".{os.sep}{self.sphinxsys_case_name} --regression=true"
        os.system(self.enter_sphinxsys_exec_folder)
        os.system(command)
        print('Simulating case is finished...')

    def run_case_with_reload(self) -> None:
        print('Start case simulation with particle reload...')
        print(self.enter_sphinxsys_exec_folder)
        command = f".{os.sep}{self.sphinxsys_case_name} --reload=true --regression=true"
        os.system(self.enter_sphinxsys_exec_folder)
        os.system(command)
        print('Simulating case is finished...')

    def read_dat_file(self):
        file = open(self.condition_file_path)
        ifconverged = file.readline(4)
        file.close()
        return ifconverged

    def clean_input_folder(folder_path, keep_file="regression_test_tool.py"):
        for filename in os.listdir(folder_path):
            file_path = os.path.join(folder_path, filename)
            if filename != keep_file and os.path.isfile(file_path):
                os.remove(file_path)
                print(f"Deleted: {file_path}")


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
        os.system(self.enter_sphinxsys_case_folder)
        os.system(command)
        print('Compiling test case is finished...')

    def test_case(self) -> None:
        print('Start test case...')
        command = "make test"
        os.system(self.enter_sphinxsys_case_folder)
        os.system(command)
        print('Testing case is finished...')

    def copy_reload(self) -> None:
        print('Start copy the reload file...')
        command = "cp -r reload bin"
        os.system(self.enter_sphinxsys_case_folder)
        os.system(command)
        print('Copying the reload file is finished...')

    def run_particle_relaxation(self) -> None:
        print('Start particle relaxation for the simulation...')
        command = f".{os.sep}{self.sphinxsys_case_name} --relax=true"
        os.chdir(self.sphinxsys_exec_path)
        os.system(command)
        print('Simulating case is finished...')

    def run_case(self) -> None:
        print('Start case simulation...')
        print(self.enter_sphinxsys_exec_folder)
        #Set the `LD_LIBRARY_PATH` to ensure it includes the path to `libsycl.so.7`.
        os.environ['LD_LIBRARY_PATH'] = '/opt/intel/oneapi/compiler/2024.0/lib:' + os.environ.get('LD_LIBRARY_PATH', '')
        command = f".{os.sep}{self.sphinxsys_case_name} --regression=true"
        os.chdir(self.sphinxsys_exec_path)
        os.system(command)
        print('Simulating case is finished...')

    def run_case_with_reload(self) -> None:
        print('Start case simulation with particle reload...')
        print(self.enter_sphinxsys_exec_folder)
        command = f".{os.sep}{self.sphinxsys_case_name} --reload=true --regression=true"
        os.chdir(self.sphinxsys_exec_path)
        os.system(command)
        print('Simulating case is finished...')

    def clean_input_folder(folder_path, keep_file="regression_test_tool.py"):
        for filename in os.listdir(folder_path):
            file_path = os.path.join(folder_path, filename)
            if filename != keep_file and os.path.isfile(file_path):
                os.remove(file_path)
                print(f"Deleted: {file_path}")

    def read_dat_file(self):
        file = open(self.condition_file_path)
        ifconverged = file.readline(4)
        file.close()
        return ifconverged
