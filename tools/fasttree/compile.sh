#!/bin/bash

# FastTree source location is version-specific (uses $BRANCH from bin/compile.sh):
#   v2.2.0 -> src/FastTree.c from the morgannprice/fasttree submodule
#   2.1.11 -> not in that repo, so build the C file embedded in the biowasm repo (../)
case "$BRANCH" in
	v2.2.0) SRC=FastTree.c ;;
	*)      SRC=../FastTree-2.1.11.c ;;
esac

emcc \
	-DNO_SSE -DUSE_DOUBLE -O3 -finline-functions -funroll-loops -Wall \
	-o ../build/fasttree.mjs "$SRC" -lm $EM_FLAGS
