#!/bin/bash

emcc wgsim.c \
    -o ../build/wgsim.mjs \
    $EM_FLAGS \
    -lm -O3 -w
