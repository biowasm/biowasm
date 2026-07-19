#!/bin/bash

# When update this list, need to update tools.json
UTILS=$(echo src/{basename,cat,comm,cut,date,df,dirname,du,echo,env,fold,head,join,ls,md5sum,paste,seq,shuf,sort,tail,tee,tr,uniq,wc}.mjs)

# Install dependencies
sudo apt-get install -y autopoint gperf help2man gettext texinfo bison
# --skip-po: avoid fetching translation (.po) files from translationproject.org (404s)
./bootstrap --skip-po
EM_GNU_NANOSLEEP

# LLVM (Emscripten 6 / clang) promotes these to hard errors by default; gnulib's older C
# (obstack.c, parse-datetime.c, ...) still trips them, so downgrade back to warnings.
WNO="-Wno-incompatible-function-pointer-types -Wno-implicit-function-declaration -Wno-implicit-int -Wno-int-conversion"

# Configure (--disable-nls to avoid Memory out of bounds error)
emconfigure ./configure \
  CC=emcc \
  --disable-nls \
  FORCE_UNSAFE_CONFIGURE=1 \
  TIME_T_32_BIT_OK=yes \
  --host=wasm32 \
  CFLAGS="-O3 $WNO"

# This program needs gcc and should not be compiled to WebAssembly
sed -i 's|$(MAKE) src/make-prime-list$(EXEEXT)|gcc src/make-prime-list.c -o src/make-prime-list$(EXEEXT) -Ilib/|' Makefile

# Pre-generate a self-contained lib/parse-datetime.c and touch it so make won't regenerate it.
# gnulib's rule runs `bison -d` then deletes parse-datetime.tab.h, but bison 3.8's -d makes the
# generated .c #include that now-deleted header. Generating without -d avoids the separate header.
( cd lib && bison parse-datetime.y \
    && sed -e 's|".*/parse-datetime.y"|"parse-datetime.y"|' parse-datetime.tab.c > parse-datetime.c \
    && rm -f parse-datetime.tab.c parse-datetime.tab.h )
touch lib/parse-datetime.c

# Make all commands and skip "man" errors
emmake make all CC=emcc -k WERROR_CFLAGS=""
emmake make $UTILS \
  CC=emcc EXEEXT=.mjs \
  CFLAGS="-O3 $EM_FLAGS $WNO" \
  -k WERROR_CFLAGS=""

# Don't throw error for unsupported features
sed -i 's/throw\("[a-z]*: TODO"\)/console.warn(\1)/g' src/*.mjs

# Move .js/.wasm files to build folder
mv $UTILS ../build/
mv ${UTILS//.mjs/.wasm} ../build/
