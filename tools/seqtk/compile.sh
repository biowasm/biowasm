#!/bin/bash

emcc src/seqtk.c \
    -o build/seqtk.html \
    $EM_FLAGS
