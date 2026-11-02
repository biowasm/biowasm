#!/bin/bash

# Remove extra space in front of Emscripten "-s" variables to avoid errors
emcmake cmake $EM_CMAKE_FLAGS -DCMAKE_CXX_FLAGS="-include alloca.h" -DCMAKE_EXE_LINKER_FLAGS="${EM_FLAGS//-s /-s} --preload-file res@/cawlign" .
emmake make cawlign
mv cawlign.{js,wasm,data} ../build/
