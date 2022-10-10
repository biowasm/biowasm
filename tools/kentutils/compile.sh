cd src/

# Compile libraries kentutils depends on (jkweb.a, jkOwnLib.a, paralib.a, embedded htslib)
emmake make topLibs

# Compile tools
cd utils/

cd bigBedToBed

emcc -O2 -std=c99 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DMACHTYPE_x86_64 -DUSE_HIC -I../../inc -I../../htslib -o bigBedToBed.o -c bigBedToBed.c

emcc -O2 -o bigBedToBed.html bigBedToBed.o ../../lib/x86_64/jkweb.a ../../htslib/libhts.a -lm \
  -s ERROR_ON_UNDEFINED_SYMBOLS=0 -s USE_ZLIB=1 --preload-file test.bb -s EXPORTED_RUNTIME_METHODS=["callMain","FS"] -s INVOKE_RUN=0 -s ALLOW_MEMORY_GROWTH=1

cd ..


cd bigWigToWig

emcc -O2 -std=c99 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DMACHTYPE_x86_64 -DUSE_HIC -I../../inc -I../../htslib -o bigWigToWig.o -c bigWigToWig.c

emcc -O2 -o bigWigToWig.html bigWigToWig.o ../../lib/x86_64/jkweb.a ../../htslib/libhts.a -lm \
  -s ERROR_ON_UNDEFINED_SYMBOLS=0 -s USE_ZLIB=1 --preload-file test.bw -s EXPORTED_RUNTIME_METHODS=["callMain","FS"] -s INVOKE_RUN=0 -s ALLOW_MEMORY_GROWTH=1

cd ..
