#!/bin/bash

EM_FLAGS="-sSTACK_SIZE=2097152 -g4 -sASSERTIONS=1"

# Remove extra space in front of Emscripten "-s" variables to avoid errors
emcmake cmake -DCMAKE_EXE_LINKER_FLAGS="${EM_FLAGS//-s /-s}" .
emmake make hyphy
mv hyphy.{js,wasm} ../build/
