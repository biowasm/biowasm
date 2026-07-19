#!/bin/bash

# Compile with SIMD
cd src/
emmake make ../../build/ssw-simd.mjs \
    CC=emcc \
    PROG=../../build/ssw-simd.mjs \
    CFLAGS="$EM_FLAGS -O3 -pipe -msse2 -msimd128"

# Compile without SIMD (use SIMDe to simulate SIMD support)
make clean
git clone https://github.com/simd-everywhere/simde.git
sed -i -e 's|include <emmintrin.h>|include "simde/simde/x86/sse2.h"|g' *.c *.h
sed -i -e 's|__m|simde__m|g' *.c *.h
sed -i -e 's|_mm_|simde_mm_|g' *.c *.h
emmake make ../../build/ssw.mjs \
    CC=emcc \
    PROG=../../build/ssw.mjs \
    CFLAGS="$EM_FLAGS -O3 -pipe"
