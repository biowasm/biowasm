#!/bin/bash

# Configure
autoheader
autoconf
emconfigure ./configure --with-htslib="../../htslib/src/" CFLAGS="-s USE_ZLIB=1 -s USE_BZIP2=1"

# Compile to WebAssembly
emmake make bcftools CC=emcc AR=emar \
    CFLAGS="-O2 -s USE_ZLIB=1 -s USE_BZIP2=1" \
    LDFLAGS="$EM_FLAGS -s ERROR_ON_UNDEFINED_SYMBOLS=0 -O2 --preload-file test/annotate.vcf@/bcftools/annotate.vcf"
