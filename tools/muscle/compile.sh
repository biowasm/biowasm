set -vexu
export EM_FLAGS=" -s USE_PTHREADS=0 $EM_FLAGS " ## disable phread
cd src/
## wget https://github.com/rcedgar/muscle/archive/refs/tags/5.1.0.tar.gz
## tar xf 5.1.0.tar.gz && cd muscle-5.1.0/src/ 
## modify omp.h + locallock.h + myutils.cpp + Makefile by myth 
emmake make  \
    CC=emcc CXX=em++

