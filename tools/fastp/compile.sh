#!/bin/bash

cp ../data/NA12878.fastq.gz testdata/
emmake make \
    TARGET="../build/fastp.js" \
    LIBS="-O3 -s USE_ZLIB=1 $EM_FLAGS --preload-file testdata@/fastp/testdata -s ASSERTIONS=1"
