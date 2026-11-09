#!/bin/bash

set -e

CMAKE_FLAGS="$EM_CMAKE_FLAGS -DCMAKE_BUILD_TYPE=Release -DHYPHY3_GIT_COMMIT=$(git rev-parse --short HEAD)"
LINKER_FLAGS="$EM_FLAGS -O3 -fexceptions"

emcmake cmake -S . -B build-wasm $CMAKE_FLAGS \
	-DHYPHY3_SIMD=OFF \
	-DCMAKE_EXE_LINKER_FLAGS="$LINKER_FLAGS"
cmake --build build-wasm --target hyphy3 -j4
cp build-wasm/hyphy3.{js,wasm} ../build/

emcmake cmake -S . -B build-wasm-simd $CMAKE_FLAGS \
	-DHYPHY3_SIMD=ON \
	-DCMAKE_EXE_LINKER_FLAGS="$LINKER_FLAGS"
cmake --build build-wasm-simd --target hyphy3 -j4
cp build-wasm-simd/hyphy3-simd.{js,wasm} ../build/
