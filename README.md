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

0. Make sure that gcc, gfrotran, wget, git, cmake are installed and uodated.

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

### Install on Windows Visual Studio

You can find a installation instruction video: https://youtu.be/m0p1nybM4v4, and install by the following steps: 
1. Install latest version Cmake, SmartGit (choose non-commercial option) binary.
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

6. Build execute able and run test cases in Visual Studio

7. Create and build your own application

        Create your own application in the cases_user in the source folder simply by copying the entire folder of a similar test case and rename and modify application files


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