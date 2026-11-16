FROM ubuntu:22.04

ARG build_with_dependencies_source=0
ARG SPH_ONLY_STATIC_BUILD=0
# 0 means auto-detect: use nproc/2 to balance build speed against OOM risk from parallel compilation
ARG build_jobs=0

ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y \
    apt-utils \
    build-essential \
    git \
    cmake \
    libgtest-dev \
    libtbb-dev \
    libboost-all-dev \
    libsimbody-dev \
    libsimbody3.6 \
    libeigen3-dev \
    libspdlog-dev \
    pybind11-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV TBB_HOME=/usr/lib/x86_64-linux-gnu
ENV BOOST_HOME=/usr/lib/x86_64-linux-gnu
ENV SIMBODY_HOME=/usr

RUN JOBS=$(( ${build_jobs} > 0 ? ${build_jobs} : $(nproc) / 2 )) && \
    mkdir -p /usr/src/googletest/build && \
    cd /usr/src/googletest/build && \
    cmake .. -DCMAKE_POSITION_INDEPENDENT_CODE=ON && \
    make -j${JOBS} && \
    make install

COPY ./ /home/SPHinXsys/
WORKDIR /home/SPHinXsys

RUN git submodule update --init
RUN rm -rf build
RUN JOBS=$(( ${build_jobs} > 0 ? ${build_jobs} : $(nproc) / 2 )) && \
    mkdir build && cd build && \
    cmake .. -DBUILD_WITH_DEPENDENCIES_SOURCE=${build_with_dependencies_source} -DSTATIC_BUILD=${SPH_ONLY_STATIC_BUILD} && \
    make -j${JOBS}