#!/bin/bash

BUILD_CONFIG=$1
if [ -z "$1" ] ; then
    BUILD_CONFIG=release
fi

BUILD_PLATFORM=$2
if [ -z "$2" ] ; then
    BUILD_PLATFORM=wasm
fi

BUILD_DIR=build_${BUILD_PLATFORM}_${BUILD_CONFIG}
INSTALL_DIR=install_${BUILD_PLATFORM}_${BUILD_CONFIG}

if [ $BUILD_PLATFORM == "wasm" ]; then

    cd ./3rd_party/wasmtbb \
        && emmake make -j$(nproc) tbb \
        && cd ../..

    echo 'Running wasm build ... '
    source /emsdk/emsdk_env.sh
    cmake -G Ninja -B$BUILD_DIR -H.\
        -DEMSCRIPTEN=1 \
        -DCMAKE_BUILD_TYPE=${BUILD_CONFIG} \
        -DCMAKE_MODULE_PATH=$EMSDK/upstream/emscripten/cmake \
        -DCMAKE_TOOLCHAIN_FILE=$EMSDK/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake 
fi

if [ $BUILD_PLATFORM == "x64" ]; then
    echo 'Running x64 build ... '
    cmake -G Ninja -B$BUILD_DIR -HSPHinXsys \
        -DBUILD_PLATFORM=${BUILD_PLATFORM} \
        -DCMAKE_BUILD_TYPE=${BUILD_CONFIG} \
        -DWASM_BUILD=0
fi

cd $BUILD_DIR && ninja && cd ..

# if [ $BUILD_PLATFORM == "wasm" ]; then
#     cp -r $BUILD_DIR/src/bin/* js/
# fi