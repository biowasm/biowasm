#!/bin/bash

emmake make CC=emcc TREE_DEST=../build/tree.js LDFLAGS="-O3 $EM_FLAGS"
