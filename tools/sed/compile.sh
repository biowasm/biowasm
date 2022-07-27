#!/bin/bash

# Install dependencies
sudo apt-get install -y pkg-config autopoint gperf help2man gettext texinfo bison
./bootstrap
EM_GNU_NANOSLEEP

# Configure
emconfigure ./configure --disable-nls

# Build
emmake make all CC=emcc -k WERROR_CFLAGS="" CFLAGS="$EM_FLAGS -O2 -s ERROR_ON_UNDEFINED_SYMBOLS=0" EXEEXT=.js
mv sed/sed.{js,wasm} ../build/
