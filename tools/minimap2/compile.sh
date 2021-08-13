cd src/

# Non-SIMD
emmake make -f Makefile.simde CFLAGS="-O2 -msimd128 -Wno-pass-failed" WASM_FLAGS="$EM_FLAGS --preload-file test/@/minimap2/" minimap2
