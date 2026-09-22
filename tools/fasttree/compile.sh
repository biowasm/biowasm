#!/bin/bash

# Compile with SIMD (FastTree requires no doubles if use SSE)
emcc \
	-O3 -finline-functions -funroll-loops -Wall \
	-msse -msimd128 \
	-o ../build/fasttree-simd.js old/FastTree-2.1.11.c -lm $EM_FLAGS

# Compile without SIMD
emcc \
	-DNO_SSE -DUSE_DOUBLE -O3 -finline-functions -funroll-loops -Wall \
	-o ../build/fasttree.js old/FastTree-2.1.11.c -lm $EM_FLAGS
