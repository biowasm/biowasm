#!/bin/bash

# Disable NLS to avoid memory out of bounds error (https://github.com/emscripten-core/emscripten/issues/11621#issuecomment-691807912)
emconfigure ./configure \
	--disable-nls \
	--disable-extensions \
	--disable-largefile


# Compile gawk to WebAssembly
cd support && emmake make && cd ..
emmake make \
	EXEEXT=.js \
	LIBS="-s ERROR_ON_UNDEFINED_SYMBOLS=0 $EM_FLAGS" \
	CFLAGS="-O2 -DNDEBUG" \
	AM_CFLAGS="" gawk.js

mv gawk.{js,wasm} ../build/
