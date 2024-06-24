FROM nvidia/cuda:12.5.0-devel-ubuntu20.04

ARG build_with_dependencies_source=0
ARG SPH_ONLY_STATIC_BUILD=0
ARG was_build=0
ARG build_with_visualization=off

ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get upgrade -y && apt-get install -y \ 
    apt-utils \
    build-essential \
    cmake \
    googletest \
    pkg-config \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/20f4e6a1-6b0b-4752-b8c1-e5eacba10e01/l_BaseKit_p_2024.0.0.49564_offline.sh
# Install Intel oneAPI
RUN sh ./l_BaseKit_p_2024.0.0.49564_offline.sh -a --silent --eula accept
# RUN source /opt/intel/oneapi/setvars.sh --include-intel-llvm

# Install the oneAPI Plugin for NVIDIA GPU
# https://developer.codeplay.com/settings/api/download/
# A new token and download URL needs to be generated. TODO: can we automate this completely?
# Remember to add version=2024.0.0
RUN wget --content-disposition 'https://developer.codeplay.com/api/v1/products/download?product=oneapi&variant=nvidia&version=2024.0.0&filters[]=linux&aat=MmZjNDU1MmNlMTVlZDE1Mzk5MGY3NTM2ZN7wAABGXejelYo4i2h1NCo3IgJJ-pf47QH1cLClaINraNLZAkfcpVTVYalyktXkInzpHrwQSpie1y23pwPkP-F5aV47'
RUN sh oneapi-for-nvidia-gpus-2024.0.0-cuda-12.0-linux.sh -y

RUN apt-get update && apt-get install -y ccache python3-dev gfortran git curl zip unzip tar \
    autoconf automake autoconf-archive

RUN cd $HOME && git clone https://www.github.com/microsoft/vcpkg.git
RUN cd $HOME && cd vcpkg && \
./bootstrap-vcpkg.sh
RUN cd $HOME && cd vcpkg && ./vcpkg install --clean-after-build \
    eigen3                          \
    tbb                             \
    boost-program-options           \
    boost-geometry                  \
    gtest                           \
    pybind11                        \
    simbody

COPY ./ /home/SPHinXsys/
WORKDIR /home/SPHinXsys
RUN rm -rf build

ENV TBB_HOME=/usr/lib/x86_64-linux-gnu
ENV BOOST_HOME=/usr/lib/x86_64-linux-gnu
ENV SIMBODY_HOME=/home/simbody
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SIMBODY_HOME/lib
ENV CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$SIMBODY_HOME/include

SHELL ["/bin/bash", "-c"]
RUN source /opt/intel/oneapi/setvars.sh --include-intel-llvm && cmake   -G "Unix Makefiles"                                                     \
    -D CMAKE_BUILD_TYPE=Release                                                 \
    -D CMAKE_C_COMPILER=icx -D CMAKE_CXX_COMPILER=icpx                          \
    -D CMAKE_TOOLCHAIN_FILE="$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake"      \
    -D CMAKE_C_COMPILER_LAUNCHER=ccache -D CMAKE_CXX_COMPILER_LAUNCHER=ccache   \
    -D SPHINXSYS_USE_SYCL=ON                                                    \
    -D SPHINXSYS_SYCL_TARGETS=nvptx64-nvidia-cuda                               \
    -S .                                                                        \
    -B ./build
RUN cmake --build build/ --target test_2d_dambreak_sycl
RUN mkdir build && cd build && cmake .. -DWASM_BUILD=${was_build} -DBUILD_WITH_DEPENDENCIES_SOURCE=${build_with_dependencies_source} -DSTATIC_BUILD=${SPH_ONLY_STATIC_BUILD} && make -j$(nproc)