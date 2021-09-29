#!/bin/bash

# Cribbed from htslib example
# we're using a development branch of htslib, so cannot use the biowasm version

sudo apt-get install -y zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev autoconf

cd src
make clean

# build htslib.a
cd htslib
make clean
autoheader
autoconf
emconfigure ./configure CFLAGS="-s USE_ZLIB=1 -s USE_BZIP2=1" --disable-lzma

# get a lot of "undefined symbol: pthread_attr_getstacksize" here
emmake make CC=emcc AR=emar \
    LDFLAGS="$EM_FLAGS -s ERROR_ON_UNDEFINED_SYMBOLS=0 -O2"
cd ..


# argp - this is a little gnarly
wget https://www.lysator.liu.se/~nisse/misc/argp-standalone-1.3.tar.gz
tar -xzvf argp-standalone-1.3.tar.gz
cd argp-standalone-1.3
wget https://raw.githubusercontent.com/Homebrew/formula-patches/b5f0ad3/argp-standalone/patch-argp-fmtstream.h
patch argp-fmtstream.h patch-argp-fmtstream.h
emconfigure ./configure
emmake make CC=emcc AR=emar \
    LDFLAGS="$EM_FLAGS -O2"
cd ..


## build program
emmake make modbam2bed CC=emcc AR=emar \
    CFLAGS="-s USE_ZLIB=1 -s USE_BZIP2=1 -fno-stack-protector -O2" \
    ARGP=argp-standalone-1.3/libargp.a \
    EXTRA_CFLAGS="-Iargp-standalone-1.3" \
    EXTRA_LDFLAGS="$EM_FLAGS --preload-file"
