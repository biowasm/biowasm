#!/bin/bash

emcc wgsim.c \
    -o ../build/wgsim.js \
    $EM_FLAGS \
    -lm -O2 -w
