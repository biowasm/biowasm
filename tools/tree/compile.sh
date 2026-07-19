#!/bin/bash

emmake make CC=emcc TREE_DEST=../build/tree.mjs LDFLAGS="-O3 $EM_FLAGS"
