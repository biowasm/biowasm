#!/bin/bash

cd src/src/

# Compile
emmake make ../../build/ssw.html \
    CC=emcc \
    PROG=../../build/ssw.html \
    CFLAGS="$EM_FLAGS -O3 -pipe -msse2 -msimd128"
