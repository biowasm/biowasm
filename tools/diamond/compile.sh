#!/bin/bash

set -e

embuilder build zlib sqlite3

CMAKE_FLAGS="$EM_CMAKE_FLAGS \
	-DCMAKE_BUILD_TYPE=Release \
	-DX86=OFF \
	-DARM=OFF \
	-DAARCH64=OFF \
	-DCROSS_COMPILE=ON \
	-DHAVE_SYSINFO=0 \
	-DSQLite3_LIBRARY= \
	-DZLIB_LIBRARY=$(em-config CACHE)/sysroot/lib/wasm32-emscripten/libz.a \
	-DZLIB_INCLUDE_DIR=$(em-config CACHE)/sysroot/include"
CXX_FLAGS="-O3 -DNDEBUG -fwasm-exceptions -sUSE_SQLITE3=1"
LINKER_FLAGS="${EM_FLAGS//-s /-s} -fwasm-exceptions -sUSE_SQLITE3=1 -sERROR_ON_UNDEFINED_SYMBOLS=0"

emcmake cmake -S . -B build-wasm $CMAKE_FLAGS \
	-DWASM_SIMD=OFF \
	-DCMAKE_CXX_FLAGS="$CXX_FLAGS" \
	-DCMAKE_EXE_LINKER_FLAGS="$LINKER_FLAGS"
emmake make -C build-wasm -j4 diamond
cp build-wasm/diamond.{js,wasm} ../build/

emcmake cmake -S . -B build-wasm-simd $CMAKE_FLAGS \
	-DWASM_SIMD=ON \
	-DCMAKE_CXX_FLAGS="$CXX_FLAGS" \
	-DCMAKE_EXE_LINKER_FLAGS="$LINKER_FLAGS"
emmake make -C build-wasm-simd -j4 diamond
cp build-wasm-simd/diamond-simd.{js,wasm} ../build/
