#!/bin/bash

# HMMER keeps Easel in a separate repo, and its build expects to find it at
# ./easel. Both are biowasm submodules (tools/hmmer/src and tools/hmmer/easel);
# copy Easel in so that `git clean -xdf` in src/ removes it again next time.
echo "Staging Easel..."
rm -rf easel
cp -r ../easel easel
rm -rf easel/.git
patch -s -p1 -d easel -i ../../patches/easel-0.49.patch

# Run autoconf to create "configure"
echo "Running autoconf..."
autoconf
( cd easel && autoconf )

# HMMER requires SIMD, so use SIMDe to simulate it
DIR_SIMDE=$(cd ../simde && pwd)
mkdir -p shim
echo "#include <simde/x86/sse.h>"  > shim/xmmintrin.h
echo "#include <simde/x86/sse2.h>" > shim/emmintrin.h
echo "#include <simde/x86/sse3.h>" > shim/pmmintrin.h
DIR_SHIM=$(pwd)/shim

# Sample data
mkdir -p test
cp tutorial/globins4.sto tutorial/globins45.fa test/

# Only included the main ones
PROGRAMS="hmmsearch hmmbuild phmmer hmmalign"

function build() {
	local suffix=$1; shift
	make clean > /dev/null 2>&1
	# ac_cv_func_popen=no: configure's probe finds popen, but Emscripten has no
	# working implementation, so HAVE_POPEN leaves hmmalign aborting on startup.
	# The only thing it costs is reading gzipped input, which needs /bin/sh anyway.
	emconfigure ./configure \
		--enable-sse --disable-threads --disable-mpi \
		ac_cv_func_popen=no \
		CC=emcc CFLAGS="-O3 $*"
	# Build only the libraries HMMER needs, not Easel's two dozen miniapps.
	emmake make -j4 -C easel libeasel.a
	emmake make -j4 -C libdivsufsort libdivsufsort.a
	emmake make -j4 -C src LDFLAGS="-s ERROR_ON_UNDEFINED_SYMBOLS=0 $EM_FLAGS --preload-file $(pwd)/test@/hmmer/"

	for program in $PROGRAMS; do
		cp src/${program}.wasm ../build/${program}${suffix}.wasm
		cp src/${program}.data ../build/${program}${suffix}.data
		sed "s|${program}\.wasm|${program}${suffix}.wasm|g; s|${program}\.data|${program}${suffix}.data|g" \
			src/${program} > ../build/${program}${suffix}.js
	done
}

echo "Compiling without SIMD"
build "" -I$DIR_SHIM -I$DIR_SIMDE -DSIMDE_ENABLE_NATIVE_ALIASES

echo "Compiling with SIMD"
build "-simd" -msimd128
