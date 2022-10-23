#!/bin/bash

emmake make release CXX=em++ LINKER=em++ LDFLAGS=
em++ -o ../build/gffread.js ../gclib/*.o gff_utils.o gffread.o -O2 $EM_FLAGS \
    --preload-file examples/transcripts.gtf@/gffread/transcripts.gtf \
    --preload-file examples/annotation.gff@/gffread/annotation.gff
