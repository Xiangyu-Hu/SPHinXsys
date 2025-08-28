FROM ubuntu:20.04

ARG build_with_dependencies_source=0
ARG SPH_ONLY_STATIC_BUILD=0
ARG was_build=0
ARG build_with_visualization=off

ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y \ 
    apt-utils \
    build-essential \
    cmake \
    googletest \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN if [ "$build_with_visualization" = on ] ; then cd /home \
    && apt-get update && apt-get install -y \ 
    libglu1-mesa-dev freeglut3-dev mesa-common-dev libxi-dev libxmu-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*; fi

RUN if [ "$build_with_dependencies_source" = 0 ] ; then cd /home \
    && apt-get update && apt-get install -y \
    libtbb-dev \
    libboost-all-dev \
    liblapack-dev \    
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && wget https://github.com/simbody/simbody/archive/Simbody-3.7.tar.gz \
    && tar xvzf Simbody-3.7.tar.gz \
    && rm Simbody-3.7.tar.gz \ 
    && mkdir /home/simbody-build && mkdir /home/simbody \
    && cd /home/simbody-build \
    && cmake /home/simbody-Simbody-3.7 -DCMAKE_INSTALL_PREFIX=/home/simbody -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBUILD_VISUALIZER=${build_with_visualization} -DBUILD_STATIC_LIBRARIES=on \
    && make -j$(nproc) \
    # && ctest -j$(nproc) \
    && make -j$(nproc) install \
    && rm -rf /home/simbody-Simbody-3.7 ; fi

ENV TBB_HOME=/usr/lib/x86_64-linux-gnu
ENV BOOST_HOME=/usr/lib/x86_64-linux-gnu
ENV SIMBODY_HOME=/home/simbody
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SIMBODY_HOME/lib
ENV CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$SIMBODY_HOME/include

COPY ./ /home/SPHinXsys/
WORKDIR /home/SPHinXsys
RUN rm -rf build
RUN mkdir build && cd build && cmake .. -DWASM_BUILD=${was_build} -DBUILD_WITH_DEPENDENCIES_SOURCE=${build_with_dependencies_source} -DSTATIC_BUILD=${SPH_ONLY_STATIC_BUILD} && make -j$(nproc)