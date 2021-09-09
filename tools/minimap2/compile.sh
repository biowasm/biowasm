cd src/

WASM_FLAGS="$EM_FLAGS --preload-file test/@/minimap2/"

# Non-SIMD
echo "Compiling without SIMD"
make clean
emmake make \
	-f Makefile.simde \
	PROGRAM="minimap2-nosimd" \
	CFLAGS="-O2 -Wno-pass-failed -Wno-return-type" \
	WASM_FLAGS="$WASM_FLAGS" \
	minimap2-nosimd

# SIMD
echo "Compiling with SIMD"
make clean
emmake make \
	-f Makefile \
	PROGRAM="minimap2" \
	sse2only=1 \
	CFLAGS="-g -O2 -msimd128" \
	WASM_FLAGS="$WASM_FLAGS -s ASSERTIONS=1" \
	minimap2

# 


exit

time minimap2 -a /minimap2/MT-human.fa /minimap2/MT-orang.fa > test.sam
time minimap2 -x map-ont -d MT-human-ont.mmi /minimap2/MT-human.fa
time minimap2 -a MT-human-ont.mmi /minimap2/MT-orang.fa > test.sam

