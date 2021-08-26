# ![](SPHINXsys/logo.png) SPHinXsys

## Description

SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle Hydrodynamics for industrial compleX systems. 
It provides C++ APIs for physical accurate simulation and aims to model coupled industrial dynamic systems including fluid, solid, multi-body dynamics 
and beyond with SPH (smoothed particle hydrodynamics), a meshless computational method using particle discretization. 
Please check the documentation of the code at https://xiangyu-hu.github.io/SPHinXsys/.
For more information on the SPHinXsys project, please check the project website: https://www.sphinxsys.org.

SPHinXsys is a multi-physics, multi-resolution SPH library. 
Although it is not a standalone application itself, 
many examples designated for the specific type of applications are provided.

## Algorithm publications

The algorithms in SPHinXsys are based on the following publications:

1. Luhui Han and Xiangyu Hu, 
"SPH modeling of fluid-structure interaction", 
Journal of Hydrodynamics, 2018: 30(1):62-69.

2. Chi Zhang and Massoud Rezavand and Xiangyu Hu, 
"Dual-criteria time stepping for weakly compressible smoothed particle hydrodynamics", 
Journal of Computational Physics 404 (2020) 109135

3. Chi Zhang et al. 
"SPHinXsys: An open-source meshless, multi-resolution and multi-physics library",
Software Impacts, 6 (2020) 100033

4. Chi Zhang, Massoud Rezavand, Xiangyu Hu,
"A multi-resolution SPH method for fluid-structure interactions",
<<<<<<< HEAD
Journal of Computational Physics,  Journal of Computational Physics, 429 (2021) 110028.
=======
Journal of Computational Physics, 2021(429) 110028
>>>>>>> e4a3571838df4a6c121d061d5da795fbc677f134

5. Chi Zhang, Yanji Wei, Frederic Dias, Xiangyu Hu,
"An efficient fully Lagrangian solver for modeling wave interaction with oscillating wave energy converter",
Ocean Engineering, 
Volume 236, (2021) 109540

6. Chi Zhang, Jianhang Wang, Massoud Rezavand, Dong Wu, Xiangyu Hu,
"An integrative smoothed particle hydrodynamics framework for modeling cardiac function",
<<<<<<< HEAD
 Computer Methods in Applied Mechanics and Engineering, 381, 113847, 2021.
=======
 Computer Methods in Applied Mechanics and Engineering, 381, 113847, 2021

 7. Yujie Zhu, Chi Zhang, Yongchuan Yu, Xiangyu Hu, 
 “A CAD-compatible body-fitted particle generator for arbitrarily complex geometry and its application to wave-structure interaction”, 
 Journal of Hydrodynamics, 33(2), 195-206, 2021.
>>>>>>> e4a3571838df4a6c121d061d5da795fbc677f134

 7. C. Zhang, M. Rezavand, Y. Zhu, Y. Yu, D. Wu, W. Zhang, J. Wang, X. Hu, "SPHinXsys: an open-source multi-physics and multi-resolution library based on smoothed particle hydrodynamics", Computer Physics Communications, 267, 108066, 2021.

 8. Yujie Zhu, Chi Zhang, Yongchuan Yu, Xiangyu Hu, "A CAD-compatible body-fitted particle generator for arbitrarily complex geometry and its application to wave-structure interaction", Journal of Hydrodynamics, 33(2), 195-206, 2021.

<<<<<<< HEAD
## Software Architecture
=======
2. Edit the CMake variables to define which dependency to use. Simbody and/or TBB can be built by the project. If one is not built by the project, install that dependency in the usual way as written before.

        Go to SPHinXsys/cmake/Dependency_free_settings.cmake
        Set BUILD_WITH_DEPENDENCIES to 1
        Set BUILD_WITH_SIMBODY to 1 if Simbody should be built by the project
        Set BUILD_WITH_ONETBB to 1 if TBB should be built by the project
        Set ONLY_3D to 1 if the 2D libraries and test cases are not needed. Note that Boost is still needed if this variable is set to 0
        Do not modify the other variables

3. Build the SPHinXsys project as described in the previous section

### Install on Windows Visual Studio

You can find a installation instruction video: https://youtu.be/m0p1nybM4v4, and install by the following steps: 
1. Install latest version Cmake, SmartGit (choose non-commercial option) binary and googel test.

        Install google test, we download the release version from the github reporsitory: <https://github.com/google/googletest/releases>, 
        build and install it. For this, you will extract the source and create a new build directory. 
        Using Cmake, you will configurate and generate a Visual Studio project. 
        Be sure that, in Cmake GUI, you have clicked the two options: build_shared_libs and install_gtest. 
        The install prefix you can choose the default one (in winodws program files and, in this case, you later need run Visual Studio as admistrator) or other new directory. 
        Open the generated project in Visual Studio, build all and install both for Debug and ReleaseWithDebugInfo targets. 
        Then, you need setup Windows system environment variables: GTEST_HOME with the vaule of the install prefix directory. 
        Also you need add the bin directory as new path. the dll files inside need to found when running the tests.

2. Build, test and install Simbody

        Downloading from https://github.com/simbody/simbody/releases
        Unpack to source folder, like: c:\simbody-source 
        Create build folder, like: c:\simbody-build
        Use Cmake, configure with option Visual Studio 2017 x64 and then Generate the solution file for VS2017 (Note that install prefix should be a file folder not in system folder. For example : C:/simbody)
        Choose RelWithDebInfo target, Right-clicking ALL_BUILD and selecting build
        Right-clicking INSTALL and selecting build.
        Choose Debug target (this important for debug with Simbody functions) 
        Right-clicking ALL_BUILD and selecting build
        Right-clicking INSTALL and selecting build.
        Set Environment Variable (User Variables) by add an entry SIMBODY_HOME to the simbody directory.
        Add the simbody\bin path to Environmental Variable (System variables)

3. Install TBB library

        Download Binary installer, actually extract the file to the assigned folder , e.g. C:/tbb_version
        Set Environment Variable (User Variables): TBB_HOME to the tbb directory
        Set the path bin\intel64\vc14 to Environmental Variable (System variables)

4. Install Boost library

        Download binary and install boost with download from https://sourceforge.net/projects/boost/files/boost-binaries/
        NOTE that you need choose right version for your visual studio. For VS 2019 you choose msvc-14.2-64, VS2017 msvc-14.1-64
        Set Environment Variable (User Variables): BOOST_HOME to its directory
        Add the Boost library (lib64 with version) path to Environmental Variable (System variables)

5. Buidling SPHinXsys project

        Clone SPHinXsys source files to local computer using SmartGit
        Remember to create a new build directory outside of the git directory to avoid upload the project files to the         
        Use Cmake to build project file
        Configure x64 build and Generate
        After configuration, one can choose debug or release mode of the project file.
        Note that you need choose the same debug or release mode in the Visual Studio as you have chosen in Cmake.
        Otherwise, it may lead to compiling issues. 

6. Build execute able and run test cases in Visual Studio

7. Create and build your own application

        Create your own application in the cases_user in the source folder simply by copying the entire folder of a similar test case and rename and modify application files

### Build with docler

Two docker files are provided:

1. Dockerfile: C++ default docker image (x86_64 ubuntu 20.04). Every libraries are installed using debian packages (Simbody is 3.6 instead of 3.7). The docker image size is smallest and suitable for docker hub or CI/CD.

command to  build: 

docker build . -t sphinxsys:latest

command to run test:

docker run sphinxsys:latest bash scripts/runTest.sh

2. dev.Dockerfile:  development packages are all installed and simbody is downloaded and compiled in the image building process. This image is too big and can not be used for github testing. Only for local development purposes. 
>>>>>>> e4a3571838df4a6c121d061d5da795fbc677f134

SPHinXsys is cross-platform can be compiled and used in Windows, Linux and Apple systems.

## Installation and tutorials

For installation, program manual and tutorials, please check https://www.sphinxsys.org/html/sphinx_index.html. 

## Contribute to SPHinXsys

You are welcomed to contribute to SPHinXsys. 
As the code is on git-hub, you can register an account there (if you do not have a github account yet) 
and fork out the SPHinXsys repository.
You can work on the forked repository and add new features, and then commit them. 

### Contribute new features

After you have finished a new feature, 
which is a commit together with a new test case which uses the new feature,
you can send a pull request to the SPHinXsys repository. 
The SPHinXsys team will review the new code and give suggestions for revision, if there is.
After the revision is done, your contribution will be merged into the master branch of the original SPHinXsys repository.
When you write the code, please remember to add your name as an author in the header files that you have contributed.

### Contribute new test cases

Test cases are very important for testing the features and showing the ability of SPHinXsys.
If you think that a cases which implemented by you are very useful for the above propose,
you can move this test into the folder cases_test in your local commit 
which is directly branched from the master branch in the original SPHinXsys repository.
After test it, you can carry out a pull request.
Then, a standard review process will be carried out before merge it to the original SPHinXsys repository.
If a new test cases is added, it will be maintained together with code by SPHinXsys team. 

### Contribute challenging benchmarks

Improving numerical algorithms is one of our main drives 
and we are happy to test challenging cases and to be informed 
if there is any other code that does better on any benchmark case.
For this, you can inform us by the contact listed below. 
If you have the case implemented, you can do just as contribute a new test case, 
except move this case to cases_user folder. 
Then, after a standard pull request and review process, 
we will merge you code to master branch and discuss with you on the issues of the case 
and work together to improve the numerical algorithm for that case. 

If you have any further question, you are also welcomed to contact xiangyu.hu@tum.de.