#!/bin/bash

# Dependencies
sudo apt-get install -y libtool

# Generate and run ./configure
autoreconf -fi
emconfigure ./configure \
    --with-oniguruma=builtin \
    --disable-maintainer-mode

# Build
emmake make EXEEXT=.mjs CFLAGS="-O3 $EM_FLAGS"
mv jq.{mjs,wasm} ../build/
