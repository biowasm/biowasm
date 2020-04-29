#!/bin/bash

cd src/

# Cleanup previous builds if any
rm bin/*
make clean
git submodule foreach 'git stash'

# Replace -lz with -s USE_ZLIB=1
find . -name 'Makefile' -exec sed -i 's/-lz/-s USE_ZLIB=1/g' {} \;
# Replace ar with $(AR) so it picks up the Emscripten AR
find . -name 'Makefile' -exec sed -i 's/ar / $(AR) /g' {} \;
# Append -s USE_ZLIB=1 in all CFLAGS (don't want to overwrite CFLAGS because some have custom "-I" statements)
find . -name 'Makefile' -exec sed -i '/^CFLAGS/ s/$/ -s USE_ZLIB=1/g' {} \;

# Build
emmake make

# Generate JS/Wasm for each tool
em++ -O3 bin/needleman_wunsch.o \
    -o ../build/needleman_wunsch.js \
    $EM_FLAGS

em++ -O3 bin/smith_waterman.o \
    -o ../build/smith_waterman.js \
    $EM_FLAGS

em++ -O3 bin/lcs.o \
    -o ../build/lcs.js \
    $EM_FLAGS
