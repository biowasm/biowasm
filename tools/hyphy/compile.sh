#!/bin/bash

# Extra compilation flags
EXTRA_FLAGS=$(cat <<EOF
    -02
    -s TOTAL_STACK=2097152
    -s EXIT_RUNTIME=0
    -s EXPORTED_RUNTIME_METHODS=['callMain','FS','PROXYFS','WORKERFS','UTF8ToString','getValue','AsciiToString']
    -fwasm-exceptions
    --preload-file res@/hyphy
EOF
)
EXTRA_FLAGS=$(echo $EXTRA_FLAGS)  # Remove whitespace

# Compile hyphy
emcmake cmake -DCMAKE_EXE_LINKER_FLAGS="$EM_FLAGS $EXTRA_FLAGS"
emmake make hyphy
mv hyphy.{js,wasm,data} ../build/
