# ![](SPHINXsys/logo.png) SPHinXsys

## Description

SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle Hydrodynamics for industrial compleX systems. 
It provides C++ APIs for physical accurate simulation and aims to model coupled industrial dynamic systems including fluid, solid, multi-body dynamics 
and beyond with SPH (smoothed particle hydrodynamics), a meshless computational method using particle discretization. 
Please check the documentation of the code at https://xiangyu-hu.github.io/SPHinXsys/.
For more information on the SPHinXsys project, please check the project website: https://www.sphinxsys.org.
For program manual and tutorials, please check https://www.sphinxsys.org/html/sphinx_index.html. 

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
Journal of Computational Physics, in press (2021)

5. Chi Zhang, Yanji Wei, Frederic Dias, Xiangyu Hu,
"An efficient fully Lagrangian solver for modeling wave interaction with oscillating wave energy converter",
arXiv:2012.05323

6. Chi Zhang, Jianhang Wang, Massoud Rezavand, Dong Wu, Xiangyu Hu,
"An integrative smoothed particle hydrodynamics framework for modeling cardiac function",
 arXiv:2009.03759

## Software Architecture

SPHinXsys is cross-platform can be compiled and used in Windows, Linux and Apple systems.

## Installation

Here, we give the instructions for installing on Ubuntu Linux, Apple OS and Windows Visual Studio.

### Installing on Ubunutu Linux and Mac OS

0. Make sure that gcc, gfrotran, wget, git, cmake and google test are installed and updated.

1. Install Boost and TBB libraries

        $ sudo apt-get install libtbb-dev
        $ sudo apt-get install libboost-all-dev
 
    and set the environment by

        $ echo 'export TBB_HOME=/usr/lib/x86_64-linux-gnu' >> ~/.bashrc
        $ echo 'export BOOST_HOME=/usr/lib/x86_64-linux-gnu' >> ~/.bashrc

2. Install Simbody library

        LAPCK library:
        $ sudo apt-get install liblapack-dev

        optinal for visualizer:
        $ sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev
        $ sudo apt-get install libxi-dev libxmu-dev
        
        Download a release version at https://github.com/simbody/simbody/releases, for example version 3.7:
        $ wget https://github.com/simbody/simbody/archive/Simbody-3.7.tar.gz
        $ tar xvzf Simbody-3.7.tar.gz
    
        Make build and install directory, and go the build folder:
        $ mkdir $HOME/simbody-build && mkdir $HOME/simbody
        $ cd $HOME/simbody-build
        
        Configure and generate Make files:
        $ cmake $HOME/simbody-Simbody-3.7 
        -DCMAKE_INSTALL_PREFIX=$HOME/simbody 
        -DCMAKE_BUILD_TYPE=RelWithDebInfo 
        -DBUILD_VISUALIZER=on (optional set to ON if simbody visualizer is going to be used) 
        -DBUILD_STATIC_LIBRARIES=on (optional, leave it off if you don't know what are you doing)

        Bulid, test and install:
        $ make -j8
        $ ctest -j8
        $ make -j8 install
    
        Allow to be found by SPHinXsys:
        Mac:
                $ echo 'export SIMBODY_HOME=$HOME/simbody' >> ~/.bash_profile
        Linux:
                $ echo 'export SIMBODY_HOME=$HOME/simbody' >> ~/.bashrc
    
        Set environment variables:
        Mac::
                $ echo 'export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$HOME/simbody-prefix/lib' >> ~/.bash_profile
        Linux:
                $ echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SIMBODY_HOME/lib' >> ~/.bashrc
                $ echo 'export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$SIMBODY_HOME/include' >> ~/.bashrc

3. Update and check environment setup before building SPHinXsys. The following commands should update the environment and report the corresponding paths.

        $ source ~/.bashrc
        $ echo $SIMBODY_HOME
        $ echo $TBB_HOME
        $ echo $BOOST_HOME

4. Build SPHinXsys

        Download the SPHinXsys from https://github.com/Xiangyu-Hu/SPHinXsys or Bitbucket if you have the link and password to the internal group repository for the newest version:
        $ git clone https://github.com/Xiangyu-Hu/SPHinXsys.git
        $ mkdir $HOME/sphinxsys-build
        $ cd $HOME/sphinxsys-build
        $ cmake /path/to/sphinxsys/source-code -DCMAKE_BUILD_TYPE=RelWithDebInfo

        You can build, test all test cases and install by
        $ make -j
        $ ctest
        
        You can play with SPHinXsys, for example run a specific test case by
        $ cd /path/to/sphinxsys-build/cases_test/test_2d_dambreak
        $ make -j 
        $ cd /bin
        $ ./test_2d_dambreak

5. Create and build your own application
  
        Create your own application in the cases_user in the source folder simply by copying the entire folder of a similar test case and rename and modify application files

        Re-run the cmake in the build folder
        $ cmake /path/to/sphinxsys/source-code -DCMAKE_BUILD_TYPE=RelWithDebInfo

        You can make and run your application
        $ cd /path/to/sphinxsys-build/cases_user/your_application_folder
        $ make -j 
        $ cd /bin
        $ ./your_application

### Installing on Ubunutu Linux using the dependency-free version

0. Note: Do not clone the submodules if you are using the default installation

1. Get all submodules, run this command in the command line of the SPHinXsys project folder

        $ git submodule update --init --recursive

2. Edit the CMake variables to define which dependency to use. Simbody and/or TBB can be built by the project. If one is not built by the project, install that dependency in the usual way as written before.

        Go to SPHinXsys/cmake/Dependency_settings.cmake
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

command to build:

docker build . -f dev.Dockerfile -t sphinxsys:dev


additional build arguement can be added to the end of the docker build command using the following syntax: --build-arg <TAG>=<VALUE>

build_with_dependencies_source=0 : default, builds with preprecompiled dependencies
build_with_dependencies_source=1 : builds dependencies together with Sphinxsys

### How to run gpuSPHinXsys cases on CUDA enabled GPUs?
The build process for GPU cases are identical to the CPU cases on all platforms, viz. Linux, Windows and Mac OSX.
However, the following notes must be considered:

        The flag -DACTIVATE_CUDA=ON must be added to let the compiler find an installed version of CUDA. 
        CUDA 11.0 or higher has to be installed.
        Choose RelWithDebInfo mode to build.
        You need to modify the "BUILD_GPU_ARCH" from the master CMake file according to the architecture of your GPU.
        Just identify the architecture of your GPU and use the proper value of BUILD_GPU_ARCH according the the available tables like this one:
        http://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/
        Two examples are 75 for Turing and 60 for Pascal.


## Contribution

Any contribution to SPHinXsys are welcome. For this, you will do the following steps:

1. Fork the project & clone locally.
2. Create an upstream remote and sync your local copy before you branch.
3. Branch for each separate piece of work.
4. Do the work, write good commit messages, and read the CONTRIBUTING file if there is one.
5. Push to your origin repository.
6. Create a new PR in GitHub.
7. Respond to any code review feedback.

For more detailed guide, please visit
https://akrabat.com/the-beginners-guide-to-contributing-to-a-github-project/