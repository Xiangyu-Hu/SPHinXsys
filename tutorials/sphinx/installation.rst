This section describes the installation procedure for the pre-compiled binaries of the SPHinXsys. 
If you want to build from source instead, please see the separate document "How to Build Simbody from Source". 
That document is included in the source zip file.

========================
How to install SPHinXsys
========================

We will start here with general information, then give platform-specific instructions in the next sections ("Linux and Mac OX S" and 
"Windows"). Be sure to check the online installation instructions that are alongside the download package for last-minute information.

Where can I find the downloads
------------------------------

Access is granted on request via email: xiangyu.hu@tum.de

What is in the download zip files
---------------------------------

The downloads include libraries, header files, documentation and example programs for SPHinXsys. 
The installation is organized as a hierarchy of directories. 
The top level directory has 5 subdirectories: lib, include, doc, and case test. 
The downloads contain static and dynamic versions of each library, in both debug and optimized form. 
All the examples are available precompiled, and source and build scripts for them are in examples/src. 

Which download do I want
-------------------------

There are separate download packages for each of the supported platforms (Linux, Mac OS X, and Windows).
There is also a source package but if you want to build from source you are reading the wrong document -- see above.

What if I have a problem
-------------------------

If you have problems, e.g., bug report and contribute to the development of SPHinXsys, 
please email to xiangyu.hu@tum.de (Xiangyu Hu) or c.zhang@tum.de (Chi Zhang).

Installation overview
---------------------

Here is the general procedure

  - Set up your machine with the required prerequisites.
  - Download the appropriate .zip package from the Downloads page.
  - Unzip into the installation directory (can be anywhere but we will suggest default locations).
  - Set path and environment variables as needed.
  - Run installation test programs to verify.
  - SimBody library 3.6.0 or higher.

The next three sections provided details specific to each of the three platforms for which we provide binaries: Linux, Mac, Windows. 
You only need to read one of these sections.


Dependencies
--------------

SPHinXsys depends on the following:

  - cross-platform building: Cmake 3.14.0 or later. See the `Cmake <https://cmake.org/>`_ webpage.
  - compiler: Visual Studio 2017 (Windows only), gcc 4.9 or later (typically on Linux), or Apple Clang (1001.0.46.3)  or later
  - google test framework
  - BOOST library (newest version)
  - TBB library (newest version)
  - Simbody library 3.6.0 or later
  - linear algebra: LAPACK 3.5.0 or later and BLAS



Installing on Unix (Linux or Mac OS X)
---------------------------------------

The only prerequisite on Mac OS X is that you have the developer kit installed, 
which you probably do already.
At a minimum, the Accelerate framework must be installed 
because that includes Lapack ad Blas libraries on which Simbody depends. 
If you download the developer kit, those libraries are installed as well.

On Linux system, LAPACK and BLAS is require, and we refer `to here
<http://www.netlib.org/lapack/>`_ and `here
<http://www.netlib.org/blas/>`_ for more details.

To install google test, in the case we have installed Cmake, if you have ROOT authority (Ubuntu):

  $ sudo apt-get install libgtest-dev
  $ cd /usr/src/gtest/
  $ sudo cmake CMakeLists.txt
  $ sudo make
  $ cd lib/
  $ sudo cp libgtest* /usr/lib/

Other wise (NO ROOT Linux):

  $ git clone https://github.com/google/googletest.git -b release-1.11.0
  $ cd googletest  
  $ mkdir build
  $ cd build
  $ cmake ../ -DCMAKE_INSTALL_PREFIX=$HOME/gtest
  $ make -j8
  $ make install

Allow to be found by cmake: 
  
  $ echo 'export GTEST_ROOT=$HOME/gtest' >> ~/.bashrc

The installation of Simbody, refers to `this link
<https://github.com/simbody/simbody#linux-or-mac-using-make>`_.
After installing Simbody correctly, set environment variable:

  -  For Mac OS X::

        $ echo 'export SIMBODY_HOME=/path/to/simbody' >> ~/.bash_profile

  -  For Linux::

		$ echo 'export SIMBODY_HOME=/path/to/simbody' >> ~/.bashrc
		$ echo 'export LIBRARY_PATH=$SIMBODY_HOME/lib64:$LIBRARY_PATH' >> ~/.bashrc
		$ echo 'export LD_LIBRARY_PATH=$LIBRARY_PATH:$LD_LIBRARY_PATH' >> ~/.bashrc
		$ echo 'export CPLUS_INCLUDE_PATH=$SIMBODY_HOME/include:$CPLUS_INCLUDE_PATH' >> ~/.bashrc

Download a release version of TBB from `their GitHub
<https://github.com/01org/tbb/releases>`_ and then unzip it to the appropriate directory on your computer and set environment variable:

  - Mac OS X::

		$ echo 'export TBB_HOME=/path/to/tbb' >> ~/.bash_profile

  - Linux::

		$ echo 'export TBB_HOME=/path/to/tbb' >> ~/.bashrc

Download a release version of BOOST from their `webpage
<https://www.boost.org/users/download/>`_ and then unzip it to the appropriate directory on your computer and set environment variable:

  - Mac OS X::

		$ echo 'export BOOST_HOME=/path/to/boost' >> ~/.bash_profile

  -  Linux::

		$ echo 'export BOOST_HOME=/path/to/boost' >> ~/.bashrc

Download the sphinxsys-linux or sphinxsys-max, and then unzip it to the appropriate directory on your computer and set environment variable \begin{itemize}

  - Mac OS X::

		$ echo 'export SPHINXSYS_HOME=/path/to/sphinxsyslibaray' >> ~/.bash_profile

  -  Linux::

		$ echo 'export SPHINXSYS_HOME=/path/to/sphinxsyslibrary' >> ~/.bashrc

and then make a build directory like sphinxsys-build with the following command:: 

    $ mkdir $HOME/sphinxsys-build
    $ cd $HOME/sphinxsys-build
    
using the following commend to build the SPHinXsys and run all the tests with the following command::

		$  cmake /path/to/sphinxsys-alpha -DCMAKE_BUILD_TYPE=RelWithDebInfo
		$ make -j
		$ ctest

You can play with SPHinXsys, for example run a specific test case by::
  
    $ cd /path/to/sphinxsys-build/cases_test/test_2d_dambreak
    $ make -j 
    $ cd /bin
    $ ./test_2d_dambreak

Right now, you can play with SPHinXsys by change the parameters. GOOD LUCK!


Installing on Ubuntu
---------------------------------------

In order for beginners to experience SPHinXsys in the Ubuntu system, 
the following installation tutorial will explain how to start from a newly 
installed Ubuntu system, install all the required programs step by step, 
and finally complete the installation of SPHinXsys.
The installation is on Ubuntu 20.04 LTS with root right.

Please note that before any installation from **apt** and **apt-get**, 
you need to run **update** or even **upgrade** command to resynchronize or update newest packages.

Press **ctrl+alt+T** on the keyboard to open the Terminal, type::

    $ sudo apt-get update
    $ sudo apt-get upgrade

In the installion process, we need somehow to use the **wget** to download 
source files from the Internet, so we need to check whether the **wget** is already in your computer by typing::

    $ wget

If you have wget installed, the system will print::

    wget: missing URL

Otherwise, it will print::

    wget command not found

Then install the **wget** on Ubuntu by typing the command below::

    $ sudo apt-get update
    $ sudo apt-get install wget

Check if **g++** is installed by typeing::

    $ g++ --version

If you didn’t install **g++**, the system will print::

    bash: g++ : command not found

Install **g++** by typeing the command::

    $ sudo apt-get install g++

Another way to install **g++** compiler is to install it as part of **build-essential** 
package. Additionally the **build-essential** package will also install additional libraries as well 
as **gcc** compiler. In most cases or if unsure this is exactly what you need::

    $ sudo apt-get install build-essential

Then check the **g++** version again, the system will print the verison of **g++**.

Make sure if you have **git** on your computer by typeing::

    $ git --version

if not, the system will print::

    bash: git: command not found

Install **git** by typing the command below::

    $ sudo apt-get install git

Then check the **git** version again.

If you would like to use debug module, check the **gdb** is in your computer or not by typeing::

    $ whereis -b gdb

Normally you will find **gdb** after you install *build-essential* package,
if not, install **gdb** by typing command below::

    $ sudo apt-get install gdb

Now we need to install **CMake**.
For a person who does not want to open the Command Line much, 
installing software present in the Ubuntu repository through the UI is very simple. 
On your Ubuntu desktop Activities toolbar, click the Ubuntu Software icon.
In the following view, click on the search icon and enter **CMake** in the search bar. 
The first package listed in the search results is the one maintained by the Snap Store. 
From the Software Manager, click on the CMake entry to **CMake** installation page and click Install button.
**CMake** will then be installed to your system.

In case of you do not find CMake in Ubuntu Software center, then install **CMake** by 
typeing those commands below in Terminal:

Install build tools and libraries that **CMake** depends on::

    $ sudo apt-get install build-essential libssl-dev

Go to the temporary directory::

    $ cd /tmp

Then, type the following command to download the source code::

    $ wget https://github.com/Kitware/CMake/releases/download/v3.20.0/cmake-3.20.0.tar.gz

Once the **tar.gz** file is downloaded, type the following command to extract it::

    $ tar -zxvf cmake-3.20.0.tar.gz

Then move to the extracted folder as follows::

    $ cd cmake-3.20.0

Finally, run the following commands to compile and install **CMake**::

    ./bootstrap

The bootstrap process may take a long time, do not interrupt it. 
When **CMake** has bootstrapped, you will get the following output:

.. figure:: figures/cmake_bootstrap_successful.png
   :width: 600 px
   :align: center

   CMake has bootstrapped

Now you can make it by typing the following command::

    $ make

And then install it as follows::

    $ sudo make install

After the **CMake** is successfully installed, you can verify its installation and 
also if the correct version is installed, through the following command::

    $ cmake --version

Now move to **LAPACK** and **BLAS**. Don't forget to move to root folder by typing::

    $ cd

Install Lapack and Blas by typing the command below::

    $ sudo apt-get install libblas-dev liblapack-dev

In the new version of **SPHinXsys**, the **Gtest** is introduced for functional test,
to intall **Gtest**, following the stpes below::

    $ sudo apt-get install libgtest-dev
    $ cd /usr/src/gtest
    $ sudo cmake CMakeLists.txt
    $ sudo make

Now you need to find where is the :code:`.a files`. Type the following command into Terminal::

    $ find . -name “libgtest*.a”

As we can see these two files were under :code:`./lib` sub folder, then type the command below::

    $ sudo cp ./lib/libgtest*.a /usr/lib

Then we make the gtest can be found by cmake::

    $ echo ‘export GTEST_ROOT=$HOME/gtest’ >> ~/.bashrc

Move to root folder. Comes to the **Boost** and **TBB** libraries::

    $ sudo apt-get install libtbb-dev
    $ sudo apt-get install libboost-all-dev

and set the environment by::

    $ echo 'export TBB_HOME=/usr/lib/x86_64-linux-gnu' >> ~/.bashrc
    $ echo 'export BOOST_HOME=/usr/lib/x86_64-linux-gnu' >> ~/.bashrc

Notice that during the installation of Boost, you might be asked to choose the aera and the city.

**SPHinXsys** use **Simbody** to calculate the multi-body dynamics, thus we need to install **Simbody**.
Here are the optional steps for visualizer of **Simbody**::

    $ sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev
    $ sudo apt-get install libxi-dev libxmu-dev

Download a release version of **Simbody** by typing the commands::

    $ wget https://github.com/simbody/simbody/archive/Simbody-3.7.tar.gz  
    $ tar xvzf Simbody-3.7.tar.gz

Make build and install directory::

    $ mkdir $HOME/simbody-build
    $ mkdir $HOME/simbody

and go the build folder::

    $ cd $HOME/simbody-build

Configure and generate Make files::

    $ cmake $HOME/simbody-Simbody-3.7 -DCMAKE_INSTALL_PREFIX=$HOME/simbody 
      -DCMAKE_BUILD_TYPE=RelWithDebInfo 
      -DBUILD_VISUALIZER=on -DBUILD_STATIC_LIBRARIES=on 

Notice that the above command is a whole command, cannot be executed separately, 
and pay attention to the space between different commands.

Then build **Simbody** by::

    $ make -j8

Note that here the :code:`-j8` means that I use 8 cores to run in parallel.
Please consider the cores on your computer to run this command.

If you want you can test **Simbody**::

    $ ctest -j8

Install **Simbdoy**::

    $ make -j8 install

Then we make **Simbody** can be found by **CMake**::

    $ echo 'export SIMBODY_HOME=$HOME/simbody' >> ~/.bashrc

Set environment variables::

    $ echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SIMBODY_HOME/lib' >> ~/.bashrc
    $ echo 'export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$SIMBODY_HOME/include' >> ~/.bashrc

If you want to use debug module of **Simbody** later in your work, 
you can modify the **DCMAKE_BUILD_TYPE** equals to **Debug**, build and install **Simbody** again by::

    $ cmake $HOME/simbody-Simbody-3.7 -DCMAKE_INSTALL_PREFIX=$HOME/simbody 
      -DCMAKE_BUILD_TYPE=Debug -DBUILD_VISUALIZER=on -DBUILD_STATIC_LIBRARIES=on 
    $ make -j8
    $ make -j8 install

Update and check environment setup before installing SPHinXsys. 
The following commands could update the environment and report the corresponding paths::

    $ source ~/.bashrc
    $ echo $SIMBODY_HOME
    $ echo $TBB_HOME
    $ echo $BOOST_HOME 

Now we can move to the last part, install **SPHinXsys**, don't forget to move to root folder.
Download the latest version of **SPHinXsys** by the command below::

    $ git clone https://github.com/Xiangyu-Hu/SPHinXsys.git

Make build directory for **SPHinXsys**::

    $ mkdir $HOME/sphinxsys-build

go to the build folder::

    $ cd $HOME/sphinxsys-build

Configure and generate Make files::

    $ cmake $HOME/SPHinXsys -DCMAKE_BUILD_TYPE=RelWithDebInfo

Notice that the path :code:`$HOME/SPHinXsys` should be path of SPHinXsys source code, you need to confirm it.

Now you can build, test all cases of **SPHinXsys** by follwoing commands::

    $ make -j7
    $ ctest
    
Please pay attention here the :code:`ctest` without parallel execution, that is becasuse the **SPHinXsys**
has the build-in function for parallel computing, if you run :code:`ctest` with :code:`-jx`, you may get some test 
cases failed.
Again, `-j7` means that I am using a 8 cores machine.  Please do not use all cores for compiling.  

or  you can choose a specific case for running, for example, the **2d_dambreak**::

    $ cd $HOME/sphinxsys-build/tests/2d_examples/test_2d_dambreak
    $ make -j8
    $ cd bin
    $ ./test_2d_dambreak

Rigth now, you have the **SPHinXsys** successfully installed in your computer, Have fun with it!

Installing on Windows
---------------------

We provide pre-built binaries for use with Visual Studio 2017. 
If you have an earlier or later version of Visual Studio, or if you are using Visual Studio Express you will likely need to build from source (not hard). See the separate build from source document referenced at the start of this chapter.

The only prerequisite on Windows is that you have a development environment (Visual Studio) and a way to unzip the .zip package. If you don’t have one already, you’ll need to install software that can perform the unzip operation. 
The installation of Simbody on Windows is refer to `Simbody's page
<https://github.com/simbody/simbody#windows-using-visual-studio>`_, 
and after that please set the system environment variable SIMBODY_HOME to the simbody prefix directory and the simbody bin path to environmental variable( System variable).


Install google test, we download the release version from the github repository: <https://github.com/google/googletest/releases>, build and install it.
For this, you will extract the source and create a new build directory. Using Cmake, you will configure and generate a Visual Studio project. 
Be sure that, in Cmake GUI,  you have clicked the two options: build_shared_libs and install_gtest. The install prefix you can choose the default one 
(in windows program files and, in this case, you later need run Visual Studio as administrator) or other new directory. 
Open the generated project in Visual Studio, build all and install both for Debug and ReleaseWithDebugInfo targets.
Then, you need setup Windows system environment variables: GTEST_HOME with the value of the install prefix directory.
Also you need add the bin directory as new path. the dll files inside need to found when running the tests.    

Install TBB, actually extract the file to the assigned folder , e.g. $C:/ tbb_2019$
set environment variable: TBB_HOME to the tbb directory, and set the path $path/to/tbb/bin/intel64/vc14$ to environmental variable (System variables).

Install boost, actually extract the file to the assigned folder, e.g. $C:/boost, and set environment: BOOST_HOME to its directory

Download the sphinxsys-win file,
and then unzip it to the appropriate directory on your computer. 
Please note you should use simple name for the  directory, 
especially not number '0', which may trigger a bug in Cmake and leads to linking error in Visual Studio. 
Set environment variable BOOST_HOME to its directory.
Using cmake for configure project as follows 


 
.. figure:: figures/cmake-sphinxsys.png
   :width: 600 px
   :align: center

   Cmake configure sphinxsys library

After configuration, one can use Visual Studio to play with SPHinXsys. GOOD LUCK!


Installing on Ubunutu Linux using the dependency-free version
-------------------------------------------------------------

Note: Do not clone the submodules if you are using the default installation!

Get all submodules, run this command in the command line of the SPHinXsys project folder::

	$ git submodule update --init --recursive

Edit the CMake variables to define which dependency to use. Simbody and/or TBB can be built by the project. 
If one is not built by the project, install that dependency in the usual way as written before.

	- Go to SPHinXsys/cmake/Dependency_settings.cmake
	- Set BUILD_WITH_DEPENDENCIES to 1
	- Set BUILD_WITH_SIMBODY to 1 if Simbody should be built by the project
	- Set BUILD_WITH_ONETBB to 1 if TBB should be built by the project
	- Set ONLY_3D to 1 if the 2D libraries and test cases are not needed. Note that Boost is still needed if this variable is set to 0
	- Do not modify the other variables

Build the SPHinXsys project as described in the previous section.

