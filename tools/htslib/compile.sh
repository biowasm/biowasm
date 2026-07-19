#!/bin/bash

# Dependencies
sudo apt-get install -yqq zlib1g-dev libbz2-dev libcurl4-gnutls-dev libssl-dev autoconf

# Compile LZMA to WebAssembly
LZMA_VERSION="5.2.5"
curl -LO "https://tukaani.org/xz/xz-${LZMA_VERSION}.tar.gz"
tar -xf xz-${LZMA_VERSION}.tar.gz
cd xz-${LZMA_VERSION}
emconfigure ./configure --disable-shared --disable-threads
emmake make -j4 CFLAGS="-Oz -fPIC -s USE_PTHREADS=0 -s EXPORT_ALL=1 -s ASSERTIONS=1"
cd -

# Run ./configure
CFLAGS="-s USE_ZLIB=1 -s USE_BZIP2=1 ${CFLAGS_LZMA}"
LDFLAGS="$LDFLAGS_LZMA"
make clean
# Use autoreconf -i so config.guess/config.sub aux files are installed (autoheader/autoconf no longer provides them).
# Use --disable-libcurl because libcurl isn't usable in wasm and would otherwise try (and fail) to compile hfile_libcurl.c.
autoreconf -i

# Always pass an explicit --host (the build machine's own triple). Newer htslib
# auto-detects it via AC_CANONICAL_HOST, but 1.10's configure.ac otherwise defaults
# host_alias to "unknown-$(uname -s)" (e.g. "unknown-Linux"), which fails config.sub.
# Passing --host is harmless for the newer versions and fixes 1.10.
emconfigure ./configure --disable-libcurl --host="$(gcc -dumpmachine)" CFLAGS="$CFLAGS" LDFLAGS="$LDFLAGS"

# Build htslib tools
TOOLS=("tabix" "htsfile" "bgzip")
for tool in ${TOOLS[@]}; do
    emmake make $tool CC=emcc AR=emar \
        CFLAGS="-O3 $CFLAGS" \
        LDFLAGS="$EM_FLAGS -O3 -s ERROR_ON_UNDEFINED_SYMBOLS=0 $LDFLAGS"
done
