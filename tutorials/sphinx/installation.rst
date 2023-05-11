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
        sudo apt install -y python3-dev

If you want a debugger for development purposes:

..  code-block:: bash

    sudo apt install -y gdb

From here, pick a workspace where the library and any dependent code will be downloaded. 
The following block will install the direct dependencies required by SPHinXsys in user-space:

..  code-block:: bash
    
    cd $HOME
    git clone https://www.github.com/microsoft/vcpkg
    cd vcpkg
    ./bootstrap-vcpkg.sh
    ./vcpkg install --clean-after-build \
        eigen3                          \
        tbb                             \
        boost-program-options           \
        boost-geometry                  \
        simbody                         \
        gtest                           \
        pybind11

Note that, some libraries can be installed simply by apt install, such as eigen3, tbb and boost.
It seems that only the simbody and gtest should be installed by using vcpkg.
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
        gtest:x64-linux-dynamic                 \
        pybind11:x64-linux-dynamic  
    cd ..

Otherwise, please refer to the official `vcpkg documentation <https://vcpkg.io/en/docs/examples/overlay-triplets-linux-dynamic.html>`_

Building SPHinXsys
^^^^^^^^^^^^^^^^^^^^^

..  code-block:: bash
    
    git clone https://github.com/Xiangyu-Hu/SPHinXsys.git sphinxsys
    cd sphinxsys
    cmake   -G "Unix Makefiles"                                                         \
            -D CMAKE_BUILD_TYPE=Release                                                 \
            -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++                           \
            -D CMAKE_TOOLCHAIN_FILE="$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake"      \
            -D CMAKE_C_COMPILER_LAUNCHER=ccache -D CMAKE_CXX_COMPILER_LAUNCHER=ccache   \
            -S .                                                                        \
            -B ./build
    cmake   --build build/ 

If you prefer to use other installed compiler in your Linux system, 
you can just change :code:`gcc` and :code:`g++` to your favorite ones. 

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
* `Python3 <https://www.python.org/>`_

Installing dependencies
^^^^^^^^^^^^^^^^^^^^^^^
Open Windows PowerShell or Git CMD, and then run the commands below one by one. 
(Before you run these commands, you can change the installation path by using the command :code:`cd ..`, etc.)

..  code-block:: pwsh
    
    git clone https://www.github.com/microsoft/vcpkg
    cd vcpkg
    .\bootstrap-vcpkg.bat
    .\vcpkg install --clean-after-build         \
        eigen3:x64-windows                      \
        tbb:x64-windows                         \
        boost-program-options:x64-windows       \
        boost-geometry:x64-windows              \
        simbody:x64-windows gtest:x64-windows
    .\vcpkg integrate install

You can also install it by using Git Bash. 
In this way, you need to change the command :code:`.\bootstrap-vcpkg.bat` to :code:`./bootstrap-vcpkg.bat` ,
i.e., you need to use the slash :code:`/` instead of the backslash:code:`\`, as follows:

..  code-block:: bash
    
    git clone https://www.github.com/microsoft/vcpkg
    cd vcpkg
    ./bootstrap-vcpkg.bat
    ./vcpkg install --clean-after-build             \
        eigen3:x64-windows                          \
        tbb:x64-windows                             \
        boost-program-options:x64-windows           \
        boost-geometry:x64-windows                  \
        simbody:x64-windows                         \
        gtest:x64-windows                           \
        pybind11:x64-windows
    ./vcpkg integrate install

Please make sure that the name of the directory for cloning vcpkg has only using plain characters, 
especially without spaces.  Otherwise, some dependent libraries, such as tbb, can not being built successfully.
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

For any other combination, please refer to the official `architecture documentation <https://vcpkg.io/en/docs/users/triplets.html>` .

If you have difficulty to install these packages, you can use the pre-compiled vcpkg files for windows as follows:

..  code-block:: pwsh
 
    git clone  https://github.com/Xiangyu-Hu/SPHinXsys_install_vcpkg_windows

To use the pre-compiled package, 
simply extract the two-volume zip file into the directory where the SPHinXsys root directory will be also located, 
then follow the rest steps to continue.

Building SPHinXsys with Visual Studio
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First clone the repository:

..  code-block:: pwsh
    
    git clone https://github.com/Xiangyu-Hu/SPHinXsys.git sphinxsys


Then, just open Visual Studio and follow the procedure given `Visual Studio document <https://learn.microsoft.com/en-us/cpp/build/cmake-projects-in-visual-studio>` .


Building SPHinXsys via cmake-gui.exe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See the figure below. Prior configuring, you must **Add Entry** and set :code:`CMAKE_TOOLCHAIN_FILE` variable 
with a :code:`FILEPATH` type pointing to :code:`<workspace>\vcpkg\scripts\buildsystems\vcpkg.cmake` .
Then, open the solution file (:code:`.sln`) generated in the :code:`build\` folder with Visual Studio.

.. figure:: figures/CMake_configure.png
   :width: 600 px
   :align: center

   CMake configures SPHinXsys library


Installing on macOS (latest) 
---------------------------------------
The procedure is given for MAC OS 13.0.1  and clang 14.0.0 (clang-1400.0.29.202).
With the assumption that you have installed Command Line Tools and python3. 

Installing dependencies
^^^^^^^^^^^^^^^^^^^^^^^

In the terminal, install the required system dependencies, homebrew, with it, 
you can install cmake, pkg-config, and others. 
Note that gfortran is essential for lapack_reference, which is needed for simbody. 

..  code-block:: bash

        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
        brew update 
        brew install cmake
        brew install pkg-config
        brew install ccache
        brew install gfortran
        brew install ninja

From here, pick a workspace where the library and any dependent code will be downloaded. 
The following block will install the direct dependencies required by SPHinXsys in user-space:

..  code-block:: bash
    
    cd $HOME
    git clone https://www.github.com/microsoft/vcpkg
    cd vcpkg
    ./bootstrap-vcpkg.sh -disableMetrics
    ./vcpkg install --clean-after-build         \
        eigen3:x64-osx                          \
        tbb:x64-osx                             \
        boost-program-options:x64-osx           \
        boost-geometry:x64-osx                  \
        simbody:x64-osx                         \
        gtest:x64-osx                           \
        pybind11:x64-osx

Building SPHinXsys
^^^^^^^^^^^^^^^^^^^^^

..  code-block:: bash
    
    git clone https://github.com/Xiangyu-Hu/SPHinXsys.git sphinxsys
    cd sphinxsys
    cmake   -G Ninja                                                                    \
            -D CMAKE_BUILD_TYPE=Release                                                 \
            -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++                     \
            -D CMAKE_TOOLCHAIN_FILE="$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake"      \
            -D CMAKE_C_COMPILER_LAUNCHER=ccache -D CMAKE_CXX_COMPILER_LAUNCHER=ccache   \
            -S .                                                                        \
            -B ./build
    cmake   --build build/ 

Running the tests and examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

They are the same as in Ubuntu Linux (See above).  

