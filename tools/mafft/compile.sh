
set -eu
export EM_FLAGS=" -s USE_PTHREADS=0 $EM_FLAGS "

cd core
emmake make CC="emcc -s TOTAL_MEMORY=360MB  -s ASSERTIONS=1 $EM_FLAGS "
cd ..
