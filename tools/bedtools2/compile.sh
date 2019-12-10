#!/bin/bash

cd src/
emmake make

# Generate .wasm/.js files
emcc obj/*.o \
    -o ../build/bedtools2.html \
    $EM_FLAGS \
    -s ERROR_ON_UNDEFINED_SYMBOLS=0
