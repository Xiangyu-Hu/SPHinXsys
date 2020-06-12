# SPHinXsys

## Description
SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle Hydrodynamics for industrial compleX systems. 
It provides C++ APIs for physical accurate simulation and aims to model coupled industrial dynamic systems including fluid, solid, multi-body dynamics 
and beyond with SPH (smoothed particle hydrodynamics), a meshless computational method using particle discretization. 
Please check the documentation of the code at https://xiangyu-hu.github.io/SPHinXsys/

SPHinXsys is a multi-physics, multi-resolution SPH library. 
Although it is not a standalone application itself, 
many examples designated for the specific type of applications are provided.

## Software Architecture
SPHinXsys is cross-platform can be compiled and used in Windows, Linux and Apple systems.

## Installation
Here, we give the instructions for installing on Ubuntu Linux, Mac OS and Windows Visual Studio.
### Installing on Ubunutu Linux and Mac OS
0. Make sure that gcc, gfrotran, wget, git, cmake are installed and uodated.

1. Install Boost and TBB libraries

	    $ sudo apt-get install libtbb-dev
	    $ sudo apt-get install libboost-all-dev
    and set the enviroment by
        
        $ echo 'export TBB_HOME=/usr/lib/x86_64-linux-gnu' >> ~/.bashrc
	    $ echo 'export BOOST_HOME=/usr/lib/x86_64-linux-gnu' >> ~/.bashrc

2. Install Simbody library

        LAPCK library:
        $ sudo apt-get install cmake liblapack-dev
        
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

3. Build SPHinXsys

        Download the SPHinXsys from https://github.com/Xiangyu-Hu/SPHinXsys or Bitbucket if you have the link and password to the internal group repository for the newest version:
		$ git clone https://github.com/Xiangyu-Hu/SPHinXsys.git
		$ mkdir $HOME/sphinxsys-build
		$ cd $HOME/sphinxsys-build
        $ cmake /path/to/sphinxsys/source-code -DCMAKE_BUILD_TYPE=RelWithDebInfo

	    Finally, you can play with SPHInXsys, for example:
		$ cd /path/to/sphinxsys-build/cases_test/test_2d_dambreak
		$ make -j 
		$ cd /bin
		$ ./test_2d_dambreak

### Install on Windows Visual Studio
1. Install lastest version Cmake, SmartGit (choose non-commercial option) binary.
2. Build, test and install Simbody

        Downloading from https://github.com/simbody/simbody/releases
        Unpack to source folder, like: c:\simobody-source 
        Create build folder, like: c:\simbody-build
        Use Cmake, configure with option Visual Studioi 2017 x64 and then Generate the solution file for VS2017 (Note that install prefix should be a file folder not in system folder. For example : C:/simbody)
        Right-clicking ALL_BUILD and selecting build
        Right-clicking INSTALL and selecting build
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
        Remember to creat a new build directory outside of the git directory to avoid upload the project files to the         
        Use Cmake to build project file
        Configure x64 build and Generate
        After configuration, one can choose debug or relaese mode of the project file.

6. Build excutable and run test cases in Visual Studio

## Contribution
Any contribution to SPHinXsys are welcome. For this, you will do the following steps:

1. Fork the project & clone locally.
2. Create an upstream remote and sync your local copy before you branch.
3. Branch for each separate piece of work.
3. Do the work, write good commit messages, and read the CONTRIBUTING file if there is one.
4. Push to your origin repository.
5. Create a new PR in GitHub.
6. Respond to any code review feedback.

For more detailed guide, please visit 
https://akrabat.com/the-beginners-guide-to-contributing-to-a-github-project/