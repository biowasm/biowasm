emmake make \
    CXX=em++ \
    EXE=../build/viral_consensus.js \
    HTSLIB_A=../../htslib/src/libhts.a \
    LIBS="-lz -lbz2" \
    CXXFLAGS="-pedantic -std=c++11 $EM_FLAGS"
