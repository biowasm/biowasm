#!/bin/bash

EM_FLAGS="-sTOTAL_STACK=2097152 -02 -sASSERTIONS=1"

# Remove extra space in front of Emscripten "-s" variables to avoid errors
emcmake cmake -DCMAKE_EXE_LINKER_FLAGS="${EM_FLAGS//-s /-s} --preload-file res@/hyphy --preload-file tests@/hyphy" .
emmake make hyphy
emmake make install
mv hyphy.{js,wasm} ../build/
