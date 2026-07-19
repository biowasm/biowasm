#!/bin/bash

# Extra compilation flags
EXTRA_FLAGS=$(cat <<EOF
    -O3
    -s TOTAL_STACK=2097152
    -s EXIT_RUNTIME=0
    -s EXPORTED_RUNTIME_METHODS=['callMain','FS','PROXYFS','WORKERFS','UTF8ToString','getValue','AsciiToString']
    -fwasm-exceptions
    --preload-file res@/hyphy
EOF
)
EXTRA_FLAGS=$(echo $EXTRA_FLAGS)  # Remove whitespace

# Compile hyphy
# clang 23 (Emscripten 6) no longer implicitly declares these: alloca.h for alloca() (avllist.cpp),
# stdlib.h for malloc()/free() (fisher_exact.cpp). Force-include both so the sources still build.
# DISABLE_FIND_PACKAGE_OpenMP: build single-threaded (no pthreads/COOP-COEP). Otherwise OpenMP
# pulls in emscripten_futex_*/__lock/__wait pthread stubs that clash at link (duplicate symbols).
# OPENMP_FOUND=FALSE: HyPhy's CMakeLists embeds if(${OPENMP_FOUND})...endif() inside
# set_target_properties(); with OpenMP disabled the var is undefined, which breaks the arg count.
# Defining it FALSE keeps the parse valid and selects the no-OpenMP branch.
emcmake cmake -DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE -DOPENMP_FOUND=FALSE -DCMAKE_CXX_FLAGS="-include alloca.h -include stdlib.h" -DCMAKE_EXE_LINKER_FLAGS="$EM_FLAGS $EXTRA_FLAGS"
emmake make hyphy
# cmake/emscripten emits the program as hyphy.js; rename to .mjs (ES6 module output).
mv hyphy.js ../build/hyphy.mjs
mv hyphy.{wasm,data} ../build/
