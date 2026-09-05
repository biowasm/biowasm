#!/bin/bash

shopt -s globstar nullglob

BUILD_FLAGS="-fwasm-exceptions -s USE_SQLITE3=1 -O3"

# configure and build native datatool binaries
./configure --with-projects=scripts/projects/datatool/project.lst --with-build-root=native
make -C native/build all_p

# configure wasm build
emconfigure ./configure AR="emar cr" \
  --with-static \
  --without-debug \
  --without-mt \
  --without-openmp \
  --without-dll \
  --without-vdb \
  --without-ngs \
  --without-gnutls \
  --without-gcrypt \
  --without-ncbicrypt \
  --without-curl \
  --without-mysql \
  --without-bdb \
  --without-sybase \
  --without-ftds \
  --without-gui \
  --without-opengl \
  --without-z \
  --without-bz2 \
  --without-lzo \
  --without-zstd \
  --without-pcre \
  --without-pcre2 \
  --without-expat \
  --without-libxml \
  --without-libxslt \
  --without-lmdb \
  --without-boost \
  --without-sablot \
  --without-xerces \
  --without-xalan \
  --without-strip \
  --host=wasm32-unknown-none \
  --with-projects=scripts/projects/blast/project.lst \
  --with-build-root=wasm

NCBI_DATATOOL_PATH="$(realpath native/bin)" make -C wasm/build CFLAGS="$BUILD_FLAGS" CXXFLAGS="$BUILD_FLAGS" LDFLAGS="$EM_FLAGS $BUILD_FLAGS -sERROR_ON_UNDEFINED_SYMBOLS=0 -sSTACK_SIZE=1048576" all_p

for f in wasm/build/app/**/*.wasm; do
    mv ${f%.wasm} ../build/$(basename $f .wasm).js
    mv $f ../build/
done
