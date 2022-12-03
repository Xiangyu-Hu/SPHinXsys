========================
How to install SPHinXsys
========================

SPHinXsys is an open-source library hosted on Github https://github.com/Xiangyu-Hu/SPHinXsys.
If you face problems about installation, want to report bugs, or any other difficulties, please report them https://github.com/Xiangyu-Hu/SPHinXsys/issues 

Requirements
------------

SPHinXsys depends on the following:

* CMake 3.16 or later
* C++17 compliant compiler

  * Visual Studio 2017 15.7 or later (Windows)
  * GCC 8 or later (Linux)
* Boost libraries
* Intel Thread Building Blocks
* Simbody 3.6.0 or later
* Eigen 3.4 or later
* GoogleTest

Installing on Ubuntu
---------------------------------------

The procedure is given for Ubuntu 20.04 LTS and considers a user having sudo privileges.
This should be identical on any more recent versions.
The home directory :code:`$HOME` is chosen as the working directory, adapt accordingly if it differs. 

Installing dependencies
^^^^^^^^^^^^^^^^^^^^^^^

In the terminal, install the required system dependencies

..  code-block:: bash

        sudo apt update
        sudo apt upgrade
        sudo apt install -y apt-utils           # package management related utility programs
        sudo apt install -y build-essential     # GCC compilation development suite and Make
        sudo apt install -y curl zip unzip tar  # when starting on a barebone Ubuntu image for bootstrapping vcpkg
        sudo apt install -y pkg-config          # for installing libraries with vcpkg
        sudo apt install -y git                 
        sudo apt install -y cmake               
        sudo apt install -y ccache              # ccache is a compiler cache. It speeds up recompilation by caching previous compilations
        sudo apt install -y python3

If you want a debugger for development purposes:

..  code-block:: bash

    sudo apt install -y gdb

From here, pick a workspace where the library and any dependent code will be downloaded. 
The following block will install the direct dependencies required by SPHinXsys in user-space:

..  code-block:: bash
    
    cd $HOME
    git clone --depth 1 --branch 2022.11.14 https://www.github.com/microsoft/vcpkg
    cd vcpkg
    ./bootstrap-vcpkg.sh
    ./vcpkg install --clean-after-build \
        eigen3                          \
        tbb                             \
        boost-program-options           \
        boost-geometry                  \
        simbody                         \
        gtest

By default, vcpkg targets the architecture *x64* and installs the *static* version of the libraries on Linux-based systems.
To install the *shared* versions, do the following:

..  code-block:: bash

    cd $HOME
    cd vcpkg
    ./vcpkg install --clean-after-build         \
        eigen3:x64-linux-dynamic                \
        tbb:x64-linux-dynamic                   \
        boost-program-options:x64-linux-dynamic \
        boost-geometry:x64-linux-dynamic        \
        simbody:x64-linux-dynamic               \
        gtest:x64-linux-dynamic
    cd ..

Otherwise, please refer to the official `documentation <https://vcpkg.io/en/docs/examples/overlay-triplets-linux-dynamic.html>`_

Building SPHinXsys
^^^^^^^^^^^^^^^^^^^^^

..  code-block:: bash
    
    git clone https://github.com/Xiangyu-Hu/SPHinXsys.git sphinxsys
    cd sphinxsys
    cmake   -G "Unix Makefiles"                                                         \
            -D CMAKE_BUILD_TYPE=Release                                                 \
            -D CMAKE_TOOLCHAIN_FILE="$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake"      \
            -D CMAKE_C_COMPILER_LAUNCHER=ccache -D CMAKE_CXX_COMPILER_LAUNCHER=ccache   \
            -S .                                                                        \
            -B ./build
    cmake   --build build/ 

Running the tests and examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To run the whole test suite:

..  code-block:: bash

    cd $HOME
    cd sphinxsys/build
    ctest -j 1 # Run each test sequentially because SPHinXsys uses all cores

    
Notice that :code:`ctest` **must run sequentially** because SPHinXsys already uses all available cores.

For running a specific case, for example, the **2d_dambreak**:

..  code-block:: bash

    cd $HOME
    cd sphinxsys/build/tests/2d_examples/test_2d_dambreak
    make -j 7 # Where 7 is the number of parallel compilation processes, adapt according to your CPU  
    cd bin
    ./test_2d_dambreak



Installing on Windows
---------------------------------------

Pre-requisites
^^^^^^^^^^^^^^^^^^^^^^^^

* Windows 7 or newer
* `Git <https://git-scm.com/download/win>`_
* `Visual Studio 2017 or newer <https://visualstudio.microsoft.com/vs/community/>`_ (mainly for `Visual Studio Build Tools <https://devblogs.microsoft.com/cppblog/updates-to-visual-studio-build-tools-license-for-c-and-cpp-open-source-projects/>`_)
* `CMake <https://cmake.org/>`_

Installing dependencies
^^^^^^^^^^^^^^^^^^^^^^^

..  code-block:: pwsh
    
    git clone --depth 1 --branch 2022.11.14 https://www.github.com/microsoft/vcpkg
    cd vcpkg
    .\bootstrap-vcpkg.bat
    .\vcpkg install --clean-after-build             \
        eigen3:x64-windows                          \
        tbb:x64-windows                             \
        boost-program-options:x64-windows           \
        boost-geometry:x64-windows                  \
        simbody:x64-windows                         \
        gtest:x64-windows
    .\vcpkg integrate install


By default, vcpkg targets the architecture *x64* and installs the *dynamic* version of the libraries on Windows system.
To install the *static* versions, replace the former install line by the following:

..  code-block:: pwsh

    .\vcpkg install --clean-after-build          \
        eigen3:x64-windows-static                \
        tbb:x64-windows-static                   \
        boost-program-options:x64-windows-static \
        boost-geometry:x64-windows-static        \
        simbody:x64-windows-static               \
        gtest:x64-windows-static

For any other combination, please refer to the official `documentation <https://vcpkg.io/en/docs/users/triplets.html>`_


Building SPHinXsys with Visual Studio
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First clone the repository:

..  code-block:: pwsh
    
    git clone https://github.com/Xiangyu-Hu/SPHinXsys.git sphinxsys


Then, just open Visual Studio and follow the procedure given `here <https://learn.microsoft.com/en-us/cpp/build/cmake-projects-in-visual-studio>`_.


Building SPHinXsys via cmake-gui.exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See the figure below. Prior configuring, you must **Add Entry** and set :code:`CMAKE_TOOLCHAIN_FILE` variable with a :code:`FILEPATH` type pointing to :code:`<workspace>\\vcpkg\\scripts\\buildsystems\\vcpkg.cmake` 
Then, open the solution file (:code:`.sln`) generated in the :code:`build\\` folder with Visual Studio.



Installing on Unix (Linux or Mac OS X)
---------------------------------------

.. warning::
    This section is **not** up-to-date. 
    It must be reworked according to the new installation procedure.


The only prerequisite on Mac OS X is that you have the developer kit installed, 
which you probably do already.
At a minimum, the Accelerate framework must be installed 
because that includes Lapack ad Blas libraries on which Simbody depends. 
If you download the developer kit, those libraries are installed as well.

On Linux system, LAPACK and BLAS is require, and we refer `to here
<http://www.netlib.org/lapack/>`_ and `here
<http://www.netlib.org/blas/>`_ for more details.

To install google test, in the case we have installed Cmake, if you have ROOT authority (Ubuntu)::

  $ sudo apt-get install libgtest-dev
  $ cd /usr/src/gtest/
  $ sudo cmake CMakeLists.txt
  $ sudo make
  $ cd lib/
  $ sudo cp libgtest* /usr/lib/

Other wise (NO ROOT Linux)::

	$ git clone https://github.com/google/googletest.git -b release-1.11.0
	$ cd googletest  
	$ mkdir build
	$ cd build
	$ cmake ../ -DCMAKE_INSTALL_PREFIX=$HOME/gtest
	$ make -j8
	$ make install

Allow to be found by cmake::

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

		$ cmake /path/to/sphinxsys-alpha -DCMAKE_BUILD_TYPE=RelWithDebInfo
		$ make -j
		$ ctest

You can play with SPHinXsys, for example run a specific test case by::
  
    $ cd /path/to/sphinxsys-build/cases_test/test_2d_dambreak
    $ make -j 
    $ cd /bin
    $ ./test_2d_dambreak

Right now, you can play with SPHinXsys by change the parameters. GOOD LUCK!
