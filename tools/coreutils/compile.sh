#!/bin/bash

# When update this list, need to update tools.json
UTILS=$(echo src/{basename,cat,comm,cut,date,df,dirname,du,echo,env,fold,head,join,ls,md5sum,paste,seq,shuf,sort,tail,tee,tr,uniq,wc}.js)

# Install dependencies
sudo apt-get install -y autopoint gperf help2man gettext texinfo bison rsync
./bootstrap
EM_GNU_NANOSLEEP

# Configure (--disable-nls to avoid Memory out of bounds error)
emconfigure ./configure \
  CC=emcc \
  --disable-nls \
  FORCE_UNSAFE_CONFIGURE=1 \
  TIME_T_32_BIT_OK=yes \
  --host=wasm32

# This program needs gcc and should not be compiled to WebAssembly
sed -i 's|$(MAKE) src/make-prime-list$(EXEEXT)|gcc src/make-prime-list.c -o src/make-prime-list$(EXEEXT) -Ilib/|' Makefile

# Bison 3.8's generated parser #includes a header that gnulib's rule deletes, so keep it in lib/
sed -i 's|&& rm -f parse-datetime.tab.c parse-datetime.tab.h|\&\& mv parse-datetime.tab.h $(top_srcdir)/lib/ \&\& rm -f parse-datetime.tab.c|' Makefile

# Make all commands and skip "man" errors
emmake make all CC=emcc -k WERROR_CFLAGS=""
emmake make $UTILS \
  CC=emcc EXEEXT=.js \
  CFLAGS="-O2 $EM_FLAGS" \
  -k WERROR_CFLAGS=""

# Don't throw error for unsupported features
sed -i 's/throw\("[a-z]*: TODO"\)/console.warn(\1)/g' src/*.js

# Move .js/.wasm files to build folder
mv $UTILS ../build/
mv ${UTILS//.js/.wasm} ../build/
