#!/bin/bash

# Remove extra space in front of Emscripten "-s" variables to avoid errors
emcmake cmake -DCMAKE_EXE_LINKER_FLAGS="${EM_FLAGS//-s /-s} --preload-file data@/hyphy" .
emmake make hyphy
mv hyphy.{js,wasm,data} ../build/
