#!/bin/bash

emcc \
	-DNO_SSE -DUSE_DOUBLE -O3 -finline-functions -funroll-loops -Wall \
	-o ../build/fasttree.js FastTree-2.1.11.c -lm $EM_FLAGS
