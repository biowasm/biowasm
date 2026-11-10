#!/bin/bash

set -e

PROGRAMS="cmalign cmbuild cmcalibrate cmemit cmfetch cmscan cmsearch cmstat cmconvert cmpress"

# Infernal's build expects Easel and HMMER in-tree, so copy them in from the repos next door
rm -rf easel hmmer simde shim
cp -r ../easel easel
cp -r ../hmmer hmmer
cp -r ../simde simde
rm -rf easel/.git hmmer/.git simde/.git
patch -s -p1 -d easel -i ../../patches/easel-infernal-1.1.5.patch
patch -s -p1 -d hmmer -i ../../patches/hmmer-infernal-1.1.5.patch

mkdir -p shim
echo "#include <simde/x86/sse.h>"  > shim/xmmintrin.h
echo "#include <simde/x86/sse2.h>" > shim/emmintrin.h
echo "#include <simde/x86/sse3.h>" > shim/pmmintrin.h

autoconf

# SIMDe lowers the SSE intrinsics to scalar C. -msimd128 instead of this is 1.6% on cmsearch
# and slower on cmcalibrate, which is not worth doubling what the CDN serves.
emconfigure ./configure \
	--host=none-none-none \
	--enable-sse --disable-threads --disable-mpi \
	ac_cv_func_popen=no \
	CC=emcc CFLAGS="-O3 -I$(pwd)/shim -I$(pwd)/simde -DSIMDE_ENABLE_NATIVE_ALIASES"
emmake make -j4 -C easel libeasel.a
emmake make -j4 -C hmmer/libdivsufsort libdivsufsort.a
emmake make -j4 -C hmmer/src libhmmer.a
emmake make -j4 -C src $PROGRAMS LDFLAGS="-s ERROR_ON_UNDEFINED_SYMBOLS=0 -s MAXIMUM_MEMORY=4GB $EM_FLAGS"
for program in $PROGRAMS; do
	cp src/$program ../build/$program.js
	cp src/$program.wasm ../build/
done
