#!/bin/bash

emcc main.c \
    -o ../build/base.js \
    $EM_FLAGS -O3
