#!/bin/bash

# Set `DYNAMIC_ZLIB` so we don't use the zlib included in the repo, which 
# otherwise causes the warning:
#   wasm-ld: warning: function signature mismatch: gzoffset
#   >>> defined as (i32) -> i32 in ./obj/fastqreader.o
#   >>> defined as (i32) -> i64 in /emsdk/upstream/emscripten/cache/sysroot/lib/wasm32-emscripten/libz.a(gzlib.c.o)
# As a result, fastp will crash if you run it on a .fastq.gz file.
cp ../data/NA12878.fastq.gz testdata/
emmake make \
    TARGET="../build/fastp.js" \
    CXXFLAGS="-std=c++11 -O3 -DDYNAMIC_ZLIB -s USE_ZLIB=1" \
    LIBS="-s USE_ZLIB=1 $EM_FLAGS --preload-file testdata@/fastp/testdata"
