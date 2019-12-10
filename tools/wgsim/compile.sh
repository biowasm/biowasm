#!/bin/bash

emcc src/wgsim.c \
    -o build/wgsim.html \
    $EM_FLAGS \
    -lm -O2 -Wall
