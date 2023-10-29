export EM_FLAGS=" -s USE_PTHREADS=0 $EM_FLAGS " ## disable phread
cd src/
emmake make  \
    CC=emcc CXX=em++

