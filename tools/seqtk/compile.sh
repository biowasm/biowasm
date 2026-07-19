#!/bin/bash

emcc seqtk.c \
    -o ../build/seqtk.mjs \
    -O3 \
    $EM_FLAGS
