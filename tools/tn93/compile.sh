#!/bin/bash

# Remove extra space in front of Emscripten "-s" variables to avoid errors
# DISABLE_FIND_PACKAGE_OpenMP: build single-threaded (no pthreads/COOP-COEP). Otherwise tn93's
# OpenMP path pulls in emscripten_futex_* / pthread runtime symbols that don't link in a
# non-pthread build.
emcmake cmake -DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE -DCMAKE_EXE_LINKER_FLAGS="${EM_FLAGS//-s /-s} --preload-file data@/tn93" .
emmake make tn93
# cmake/emscripten emits the program as tn93.js; rename to .mjs (ES6 module output).
mv tn93.js ../build/tn93.mjs
mv tn93.{wasm,data} ../build/
