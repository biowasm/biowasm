#!/bin/bash

cd src/
emcc seqtk.c \
    -o ../build/seqtk.html \
    $EM_FLAGS
