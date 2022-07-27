#!/bin/bash

# Install dependencies
sudo apt-get install -y pkg-config autopoint gperf help2man gettext texinfo bison
./bootstrap
EM_GNU_NANOSLEEP

# Avoid grep --help showing `Usage: (null) ...`
sed -i 's/getprogname ()/"grep"/g' src/grep.c

# Configure
emconfigure ./configure --disable-nls

# Run Makefile on each subfolder. Don't error on undefined _splice function that isn't used
emmake make all CC=emcc -k WERROR_CFLAGS="" CFLAGS="$EM_FLAGS -O2 -s ERROR_ON_UNDEFINED_SYMBOLS=0" EXEEXT=.js
mv src/grep.{js,wasm} ../build/
