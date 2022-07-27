#!/bin/bash

emcc seqtk.c \
    -o ../build/seqtk.js \
    -O2 \
    $EM_FLAGS
