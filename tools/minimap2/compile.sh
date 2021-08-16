cd src/

WASM_FLAGS="$EM_FLAGS --preload-file test/@/minimap2/"

# # Non-SIMD
# echo "Compiling without SIMD"
# make clean
# emmake make \
# 	-f Makefile.simde \
# 	PROGRAM="minimap2-nosimd" \
# 	CFLAGS="-O2 -msimd128 -Wno-pass-failed" \
# 	WASM_FLAGS="$WASM_FLAGS" \
# 	minimap2-nosimd

# SIMD
echo "Compiling with SIMD"
make clean
emmake make \
	-f Makefile \
	PROGRAM="minimap2" \
	sse2only=1 \
	CFLAGS="-O2 -msimd128" \
	WASM_FLAGS="$WASM_FLAGS -s ERROR_ON_UNDEFINED_SYMBOLS=0" \
	minimap2

exit

minimap2 -v 100 -a /minimap2/MT-human.fa /minimap2/MT-orang.fa > test.sam
