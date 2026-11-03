#!/bin/bash

set -e

# Match upstream's default builds: 2.1.11 is single precision (SSE3), 2.2.0 is double precision
if [[ "$BRANCH" == "main" ]]; then
	emcc -DNO_SSE -O3 -finline-functions -funroll-loops -Wall \
		-o ../build/fasttree.js old/FastTree-2.1.11.c -lm $EM_FLAGS $EM_FLAGS_SIMD
	emcc -msse3 -msimd128 -O3 -finline-functions -funroll-loops -Wall \
		-o ../build/fasttree-simd.js old/FastTree-2.1.11.c -lm $EM_FLAGS
else
	emcc -O3 -finline-functions -funroll-loops -Wall \
		-o ../build/fasttree.js FastTree.c -lm $EM_FLAGS $EM_FLAGS_SIMD
	emcc -msimd128 -O3 -finline-functions -funroll-loops -Wall \
		-o ../build/fasttree-simd.js FastTree.c -lm $EM_FLAGS
fi
