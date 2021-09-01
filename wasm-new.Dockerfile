
FROM emscripten/emsdk:latest

ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y \ 
    apt-utils \
    cmake \
    ninja-build \
    ccache \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV CCACHE_DIR /app/.ccache

WORKDIR /app