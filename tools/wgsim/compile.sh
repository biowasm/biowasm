#!/bin/bash

cd src/
emcc wgsim.c \
    -o ../build/wgsim.html \
    $EM_FLAGS \
    -lm -O2 -Wall
