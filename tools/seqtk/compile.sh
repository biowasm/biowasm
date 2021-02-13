#!/bin/bash

cd src/
emcc seqtk.c \
    -o ../build/seqtk.html \
    -O2 \
    $EM_FLAGS
