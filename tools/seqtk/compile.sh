#!/bin/bash

emcc seqtk.c \
    -o ../build/seqtk.mjs \
    -O2 \
    $EM_FLAGS
