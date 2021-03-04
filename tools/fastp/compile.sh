#!/bin/bash

cd src/
curl -o "testdata/NA12878.fastq.gz" "https://sandbox.bio/data/NA12878.30k.fastq.gz"
emmake make TARGET="fastp.o" LIBS="-s USE_ZLIB=1"
em++ -O3 fastp.o \
    -o ../build/fastp.html \
    $EM_FLAGS \
    --preload-file testdata@/fastp/testdata \
    -s ASSERTIONS=1
