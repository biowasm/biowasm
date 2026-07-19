#!/bin/bash

# Install dependencies
sudo apt-get install -y pkg-config autopoint gperf help2man gettext texinfo bison
# --skip-po: avoid fetching translation (.po) files from translationproject.org (404s)
./bootstrap --skip-po
EM_GNU_NANOSLEEP

# LLVM (Emscripten 6 / clang) promotes these to hard errors by default; gnulib's older C
# still trips them, so downgrade back to warnings.
WNO="-Wno-incompatible-function-pointer-types -Wno-implicit-function-declaration -Wno-implicit-int -Wno-int-conversion"

# Configure
emconfigure ./configure --disable-nls CFLAGS="-O3 $WNO"

# Build
emmake make all CC=emcc -k WERROR_CFLAGS="" CFLAGS="$EM_FLAGS -O3 -s ERROR_ON_UNDEFINED_SYMBOLS=0 $WNO" EXEEXT=.mjs
mv sed/sed.{mjs,wasm} ../build/
