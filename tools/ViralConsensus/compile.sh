version=$(git describe --tags)
is_less_than_version_0_0_4=$([[ "0.0.4" != $(echo -e "$version\n0.0.4" | sort -V | head -n 1) ]] && echo true)

# ViralConsensus < 0.0.4
if [[ "$is_less_than_version_0_0_4" == "true" ]]; then
    # Because Emscripten does not support the -fpermissive flag, we need to modify htslib source code so that it compiles without errors.
    sed -i 's/unsigned char \*tmp = realloc(b->data, len);/unsigned char \*tmp = static_cast<unsigned char \*>(realloc(b->data, len));/g' htslib/cram/cram_io.h

    emmake make \
        CXX=em++ \
        EXE=../build/viral_consensus.js \
        HTSLIB_A=../../htslib/src/libhts.a \
        LIBS="-lz -lbz2 -llzma $LDFLAGS_LZMA" \
        CXXFLAGS="-std=c++11 -fpermissive $EM_FLAGS $CFLAGS_LZMA"

# ViralConsensus >= 0.0.4 doesn't rely on local htslib
else 
    emmake make \
        CXX=em++ \
        EXE=../build/viral_consensus.js \
        INCLUDE="-I../../htslib/src/" \
        LIBS="-O2 -lhts -L../../htslib/src/ -lz -lbz2 -llzma $LDFLAGS_LZMA" \
        CXXFLAGS="-Wall -pedantic -std=c++11 $EM_FLAGS -s ERROR_ON_UNDEFINED_SYMBOLS=0"
fi
