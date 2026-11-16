====================================
Set regression test for a new case
====================================

If you have built a new test case, and you want to put it into the source of the SPHinXsys. 
You need to generate a converged database for regression testing first. 


Observation variables
-------------------------------

Set up probes in some typical locations to monitor representative variables, 
such as the displacement of the beam tip, velocity profile, total force, maximum speed, etc. 
You can choose more than one such variable to set monitoring points.


Change the monitoring class 
------------------------------

Change the monitoring class :code:`ObservedQuantityRecording`
and/or :code:`BodyReducedQuantityRecording`  to be inherited by regression testing class 
:code:`RegressionTestTimeAveraged` or :code:`RegressionTestEnsembleAveraged`
or :code:`RegressionTestDynamicTimeWarping`.
The curve type of the monitored time-series data determines which method is more suitable to use.


Generate database 
------------------------------

1) After the main loop, call the method :code:`generateDataBase` of the object defined by 
the above-mentioned three regression test classes with the threshold value, like the following formats. 
The time-averaged method and ensemble-averaged method need two threshold values for both mean and variance, 
and its dimension should be the same as the variables. 
The dynamic time warping method just needs one threshold value no matter the dimension.

:code:`write_total_force_on_inserted_body.generateDataBase({ 0.005, 0.005 }, { 0.005, 0.005 });`

:code:`write_beam_tip_displacement.generateDataBase({ 0.01, 0.05 }, { 0.01, 0.05 });`

:code:`write_recorded_water_pressure.generateDataBase(0.005);`


2) Make a python script referring to the :code:`test_2d_dambreak` test case, 
and add the following command into the :code:`CMakeLists.txt`, 
which can copy the script into the input folder after making file with the CMake. 
The name of the script is generally :code:`regression_test_tool.py`.
:code:`file(MAKE_DIRECTORY ${BUILD_INPUT_PATH})`

:code:`execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_INPUT_PATH})`

:code:`file(COPY${CMAKE_CURRENT_SOURCE_DIR}/regression_test_tool/regression_test_tool.py DESTINATION ${BUILD_INPUT_PATH})`


3) Find file regression_test_tool.py 
in the path :code:`bin/input`, and open terminal, 
then type :code:`python3 regression_test_tool.py` (ONLY FOR Linux). 
It will automatic run multiply times to generate the converged database.
Generally, the mean, variance, or distance can be converged within 100 run times. 
If not, you can adjust thresholds slightly. 
The execution will stop at 200 times if it still not converged.


4) Finally, you will get the converged database with the XML format and put them into the source like the fourth step,
and add the commend to copy database files to :code:`CMakeLists.txt` for your case for building the project of the case.
You can refer the existing  :code:`CMakeLists.txt` for test case with regression test.
Then those databases can be used for regression testing. Please continue to read.


Test with the database 
------------------------------

After the main loop, call the method :code:`newResultTest` of the object defined by the regression test classes, 
like the following formats.

:code:`write_total_force_on_inserted_body.newResultTest();`

:code:`write_beam_tip_displacement.newResultTest();` 

:code:`write_recorded_water_pressure.newResultTest();` 


Other comments:
------------------------------

1)	Don’t output the monitoring data of the 0th iteration before and during the main loop at the same time. 
It will crash the XML file.

2)	The threshold value is generally set to be 0.005, but it can also be adjusted accordingly.

3)	If you are failed to pass all the cases, don’t be worry, and please try it again. 
Some cases are not stable for computing results.

4)	When you build the regression testing environment, please don’t put the code into the “src” folder, 
and just put it in the case folder, referring to the existed test cases.

5)	Add the following command to the :code:`CMakeList.txt`, 
and the resulting path of the testing will be the same as the execution path.

:code:`set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")`

6)	For the test cases that need to do the particle relaxation, 
their regression testing tool may have some differences. 
Referring to the :code:`test_2d_fsi`, whose script has the following two lines to first generate the relaxed particles.

:code:`sphinxsys.test_case()`

:code:`sphinxsys.copy_reload()`

7)	When you add a test to the CMakeList.txt, please add a working directory to the test. 

:code:`add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})`

