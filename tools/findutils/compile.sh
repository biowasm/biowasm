#!/bin/bash

# Install dependencies
sudo apt-get install -y pkg-config autopoint gperf help2man gettext texinfo bison

# Configure
./bootstrap
EM_GNU_NANOSLEEP
EM_GNU_STRCASESTR_LINEAR
emconfigure ./configure --disable-nls --disable-year2038

# Build
emmake make all CC=emcc -k WERROR_CFLAGS="" CFLAGS="$EM_FLAGS -O2 -s ERROR_ON_UNDEFINED_SYMBOLS=0" EXEEXT=.js
mv find/find.{js,wasm} ../build/
