#!/bin/bash

# Remove extra space in front of Emscripten "-s" variables to avoid errors
# -include alloca.h: newer libc/clang no longer declare alloca() implicitly (alignment.cpp uses it)
# DISABLE_FIND_PACKAGE_OpenMP: build single-threaded (no pthreads/COOP-COEP), else OpenMP pulls in
# emscripten_futex_* / pthread runtime symbols that don't link in a non-pthread build.
emcmake cmake -DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE -DCMAKE_CXX_FLAGS="-include alloca.h" -DCMAKE_EXE_LINKER_FLAGS="${EM_FLAGS//-s /-s} --preload-file res@/cawlign" .
emmake make cawlign
mv cawlign.{js,wasm,data} ../build/
