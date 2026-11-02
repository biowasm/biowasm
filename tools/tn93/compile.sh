#!/bin/bash

# Remove extra space in front of Emscripten "-s" variables to avoid errors
emcmake cmake $EM_CMAKE_FLAGS -DCMAKE_CXX_FLAGS="-Wno-error" -DCMAKE_EXE_LINKER_FLAGS="${EM_FLAGS//-s /-s} --preload-file data@/tn93" .
emmake make tn93
mv tn93.{js,wasm,data} ../build/
