#!/bin/bash

# Configure
./autogen.sh
emconfigure ./configure \
    --with-hts=../../../htslib/src \
    CPPFLAGS="-s USE_ZLIB=1 -s USE_BZIP2=1"

# Need `-s ERROR_ON_UNDEFINED_SYMBOLS=0` to avoid error `undefined symbol: BZ2_bzBuffToBuffCompress (referenced by top-level compiled C/C++ code)`
# -llzma: link htslib's liblzma (from the xz build) explicitly — the autotools link line doesn't
# add it, so libhts.a's lzma_* symbols would otherwise be undefined (fatal under Emscripten 6's -Werror).
# ivar's autotools build links from a nested dir, so LDFLAGS_LZMA's relative -L resolves wrong;
# use an absolute path to liblzma (xz is already built by the htslib dependency at this point).
LZMA_LIBDIR_ABS="$(realpath "${LDFLAGS_LZMA#-L}")"
emmake make \
    EXEEXT=.mjs \
    CFLAGS="-O3 -s USE_ZLIB=1 -s USE_BZIP2=1 $CFLAGS_LZMA" \
    LDFLAGS="$EM_FLAGS -s ERROR_ON_UNDEFINED_SYMBOLS=0 -s USE_BZIP2=1 -O3 -L${LZMA_LIBDIR_ABS} -llzma"

mv src/ivar.{mjs,wasm} ../build/
