#!/bin/bash

./autogen.sh
emconfigure ./configure
emmake make EXEEXT=".mjs" CFLAGS="-O3 -s ERROR_ON_UNDEFINED_SYMBOLS=0 $EM_FLAGS"
cp src/kalign.* ../build/
