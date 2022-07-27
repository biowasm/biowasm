#!/bin/bash

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
