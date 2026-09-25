#!/bin/bash

LINKER_FLAGS="${EM_FLAGS//-s /-s} --preload-file $PWD/data@/parasail/"

# Compile with SIMD: build SSE2 kernels and lower them to WebAssembly SIMD.
# SSE2_C_FLAGS/HAVE_SSE2 are set explicitly because cmake/FindSSE2.cmake probes
# for an x86 compiler flag that emcc does not advertise.
echo "Compiling with SIMD"
rm -rf CMakeCache.txt CMakeFiles
emcmake cmake . \
	-DCMAKE_BUILD_TYPE=Release \
	-DHAVE_SSE2=1 \
	-DSSE2_C_FLAGS="-msse2 -msimd128" \
	-DCMAKE_C_FLAGS="-O3 -msimd128" \
	-DCMAKE_CXX_FLAGS="-O3 -msimd128" \
	-DPARASAIL_ALIGNER_NAME="parasail_aligner-simd" \
	-DCMAKE_EXE_LINKER_FLAGS="$LINKER_FLAGS"
emmake make parasail_aligner
mv parasail_aligner-simd.{js,wasm,data} ../build/

# Compile without SIMD: with no vector implementation available, parasail falls
# back to its own serial reference implementations. The output name has to be set
# at build time because Emscripten bakes the .data filename into the glue code.
echo "Compiling without SIMD"
rm -rf CMakeCache.txt CMakeFiles
emcmake cmake . \
	-DCMAKE_BUILD_TYPE=Release \
	-DHAVE_SSE2=0 \
	-DSSE2_C_FLAGS="" \
	-DCMAKE_C_FLAGS="-O3" \
	-DCMAKE_CXX_FLAGS="-O3" \
	-DCMAKE_EXE_LINKER_FLAGS="$LINKER_FLAGS"
emmake make parasail_aligner
mv parasail_aligner.{js,wasm,data} ../build/
