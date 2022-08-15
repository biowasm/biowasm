#!/bin/bash

# Dependencies
sudo apt-get install -y libtool

# Generate and run ./configure
autoreconf -fi
emconfigure ./configure \
    --with-oniguruma=builtin \
    --disable-maintainer-mode

# Build
emmake make EXEEXT=.js CFLAGS="-O2 $EM_FLAGS"
mv jq.{js,wasm} ../build/
