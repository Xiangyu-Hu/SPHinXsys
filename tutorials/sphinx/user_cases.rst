===================================
Creating your own cases
===================================

There are two ways to use the library for your own problems.
One is to install the library to a destination folder and use it as a third part library.
The second and easiest way is creating new user examples within the project.

This can be done by simply copy a test example (with the entire folder, such as test_2d_dambreak) 
into the folder called user_examples. Please copy a 2D example if your own cases is 2D, 
a 3D example if your own case is 3D too. 

Then, you can modify the names of the folder and cpp file to create a new case. 
After that, you can modify the cpp or header files for your own application. 
Finally, you rerun Cmake for the entire project so that your own case will be included.

After that you can build your own case by the following command in the terminal:

..  code-block:: bash
    
    cmake   --build build/ --target user_example_your_case_name

Please replace your_case_name by the name of your own case folder.