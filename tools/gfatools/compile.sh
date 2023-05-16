emmake make \
    CC=emcc \
    CFLAGS="-std=c99 -O2 -s USE_ZLIB=1" \
    LIBS="$EM_FLAGS --preload-file test@/gfatools/test"
