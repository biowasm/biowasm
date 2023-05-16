sudo apt-get install -y zlib1g-dev

emmake make \
    CC=emcc \
    LIBS="$EM_FLAGS --preload-file test@/gfatools/test"
