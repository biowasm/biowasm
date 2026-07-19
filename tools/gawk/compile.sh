#!/bin/bash

# gawk ships pre-generated autotools files. `git checkout` skews their mtimes so the build
# thinks they're stale and re-runs automake, which then reports a version mismatch. Touch the
# generated files so they're newer than their sources and regeneration is skipped.
touch aclocal.m4 configure config.h.in
find . -name "Makefile.in" -exec touch {} +

# Disable NLS to avoid memory out of bounds error (https://github.com/emscripten-core/emscripten/issues/11621#issuecomment-691807912)
emconfigure ./configure \
	--disable-nls \
	--disable-extensions \
	--disable-largefile


# Compile gawk to WebAssembly
cd support && emmake make && cd ..
emmake make \
	EXEEXT=.mjs \
	LIBS="-s ERROR_ON_UNDEFINED_SYMBOLS=0 $EM_FLAGS" \
	CFLAGS="-O3 -DNDEBUG" \
	AM_CFLAGS="" gawk.mjs

mv gawk.{mjs,wasm} ../build/
