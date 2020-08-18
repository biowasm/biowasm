#!/bin/bash

cd src/

# Remove large files (we'll pre-load the rest of the files as examples)
rm -rf ./test/intersect/sortAndNaming/bigTests

# Generate obj/*.o files
make clean
emmake make

# Generate .wasm/.js files
emcc obj/*.o \
    -o ../build/bedtools2.html \
    -O2 \
    --preload-file test@/bedtools2/test \
    $EM_FLAGS \
    -s ERROR_ON_UNDEFINED_SYMBOLS=0
