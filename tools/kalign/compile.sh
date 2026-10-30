#!/bin/bash

# autogen.sh won't generate ./configure without libtoolize
sudo apt-get install -y libtool

./autogen.sh

# Emscripten now ships libopenmp, so configure finds it and kalign aborts spawning threads
emconfigure ./configure --without-openmp

emmake make EXEEXT=".js" CFLAGS="-O3 -s ERROR_ON_UNDEFINED_SYMBOLS=0 $EM_FLAGS"
cp src/kalign.* ../build/
