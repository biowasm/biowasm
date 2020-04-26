#!/bin/bash

cd src/
emmake make TARGET="fastp.o" LIBS="-s USE_ZLIB=1"
em++ -O3 fastp.o \
    -o ../build/fastp.js \
    $EM_FLAGS \
    --preload-file testdata@/fastp/testdata \
    -s ASSERTIONS=1
