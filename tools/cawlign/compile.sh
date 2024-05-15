#!/bin/bash

# Remove extra space in front of Emscripten "-s" variables to avoid errors
emcmake cmake -DCMAKE_EXE_LINKER_FLAGS="${EM_FLAGS//-s /-s} --preload-file test@/cawlign --preload-file res@/cawlign" .
emmake make cawlign
mv cawlign.{js,wasm,data} ../build/
