set -eu
mkdir -p testdata
cp ../data/* testdata/
em++ -static -O3 -ffast-math -lm $EM_FLAGS --preload-file  testdata@/TMalign/testdata  -o ../build/TMalign.js TMalign.cpp  

