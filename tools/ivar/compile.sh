#!/bin/bash

# Configure
./autogen.sh
emconfigure ./configure \
    --with-hts=../../../htslib/src \
    CPPFLAGS="-s USE_ZLIB=1 -s USE_BZIP2=1"

# Need `-s ERROR_ON_UNDEFINED_SYMBOLS=0` to avoid error `undefined symbol: BZ2_bzBuffToBuffCompress (referenced by top-level compiled C/C++ code)`
emmake make \
    EXEEXT=.js \
    CFLAGS="-O2 -s USE_ZLIB=1 -s USE_BZIP2=1" \
    LDFLAGS="$EM_FLAGS -s ERROR_ON_UNDEFINED_SYMBOLS=0 -O2"

mv src/ivar.{js,wasm} ../build/
