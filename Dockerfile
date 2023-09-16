FROM ubuntu:20.04

ARG build_with_dependencies_source=0
ARG SPH_ONLY_STATIC_BUILD=0

ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y \ 
    apt-utils \
    bash \
    build-essential \
    cmake \
    libgtest-dev \
    libtbb-dev \
    libboost-all-dev \
    libsimbody-dev \
    libsimbody3.6 \
    libeigen3-dev \
    python3-pybind11 \
    git \
    sudo \   
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV TBB_HOME=/usr/lib/x86_64-linux-gnu
ENV BOOST_HOME=/usr/lib/x86_64-linux-gnu
ENV SIMBODY_HOME=/usr

COPY ./ /home/SPHinXsys/
WORKDIR /home/SPHinXsys

RUN cd /usr/src/googletest && \
    sudo cmake . && \
    sudo cmake --build . --target install && \
    cd /home/SPHinXsys

RUN git submodule update --init
RUN rm -rf build
RUN mkdir build && cd build && cmake .. -DBUILD_WITH_DEPENDENCIES_SOURCE=${build_with_dependencies_source} -DSTATIC_BUILD=${SPH_ONLY_STATIC_BUILD}

ENTRYPOINT [ "/bin/bash" ]