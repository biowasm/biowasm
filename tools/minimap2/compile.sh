#!/bin/bash

WASM_FLAGS="$EM_FLAGS --preload-file test/@/minimap2/"

# Non-SIMD
echo "Compiling without SIMD"
make clean
emmake make \
	-f Makefile.simde \
	PROGRAM="minimap2" \
	CFLAGS="-O2 -Wno-pass-failed -Wno-return-type" \
	WASM_FLAGS="$WASM_FLAGS" \
	minimap2

# SIMD
echo "Compiling with SIMD"
make clean
emmake make \
	-f Makefile \
	PROGRAM="minimap2-simd" \
	CFLAGS="-O2 -Wno-return-type -Wno-unused-command-line-argument -msimd128" \
	WASM_FLAGS="$WASM_FLAGS" \
	minimap2-simd
