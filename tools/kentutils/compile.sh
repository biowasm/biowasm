# The link step below looks for lib/x86_64, so build the libs there regardless of host arch
export MACHTYPE=x86_64

# Install dependencies
sudo apt-get install -y libpng-dev libmysqlclient-dev

# Compile libraries kentutils depends on (jkweb.a, jkOwnLib.a, paralib.a, embedded htslib)
cd src/
emmake make topLibs

# Compile tools
cd utils/
TOOLS=("bigBedToBed" "bigBedInfo" "bigWigToWig" "bigWigInfo")
for TOOL in ${TOOLS[@]};
do
  cd $TOOL

  emcc -O2 -o $TOOL.o -c $TOOL.c \
    -std=c99 -I../../inc -I../../htslib \
    -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DMACHTYPE_x86_64

  emcc -O2 -o ../../../../build/$TOOL.js $TOOL.o \
    ../../lib/x86_64/jkweb.a ../../htslib/libhts.a -lm \
    -s ERROR_ON_UNDEFINED_SYMBOLS=0 $EM_FLAGS

  cd ..
done

