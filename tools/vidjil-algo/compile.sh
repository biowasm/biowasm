#!/usr/bin/env bash

EM_FLAGS="$EM_FLAGS --preload-file ../demo/Demo-X5.fa@demo.fa"

make clean
make -C algo clean
make germline
make demo
emmake make -C algo EM_FLAGS="$EM_FLAGS" wasm 

mv vidjil-algo.* ../build
