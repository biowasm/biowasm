# Because Emscripten does not support the -fpermissive flag, we need to modify htslib source code so that it compiles without errors.
sed -i 's/unsigned char \*tmp = realloc(b->data, len);/unsigned char \*tmp = static_cast<unsigned char \*>(realloc(b->data, len));/g' htslib/cram/cram_io.h

emmake make \
    CXX=em++ \
    EXE=../build/viral_consensus.js \
    HTSLIB_A=../../htslib/src/libhts.a \
    LIBS="-lz -lbz2" \
    CXXFLAGS="-std=c++11 -fpermissive $EM_FLAGS"
