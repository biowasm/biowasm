#!/bin/bash

cd src/
curl -o "testdata/NA12878.fastq.gz" "https://sandbox.bio/data/NA12878.30k.fastq.gz"
emmake make \
    TARGET="../build/fastp.html" \
    LIBS="-s USE_ZLIB=1 $EM_FLAGS --preload-file testdata@/fastp/testdata -s ASSERTIONS=1"
