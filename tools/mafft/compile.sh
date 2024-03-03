#export EM_FLAGS='-s USE_ZLIB=1 -s INVOKE_RUN=0 -s FORCE_FILESYSTEM=1 -s EXPORTED_RUNTIME_METHODS=["callMain","FS","PROXYFS","WORKERFS"] -s MODULARIZE=1 -s ENVIRONMENT="web,worker" -s ALLOW_MEMORY_GROWTH=1 -lworkerfs.js -lproxyfs.js' ## add -s ERROR_ON_UNDEFINED_SYMBOLS=0 to skip error:  undefined symbol:  参考是的是https://groups.google.com/g/emscripten-discuss/c/HSRgQiIq1gI/m/Kt9oFWHiAwAJ

set -eu
export EM_FLAGS=" -s USE_PTHREADS=0 $EM_FLAGS "
cd ../
if [ ! -d "mafft-v7.520" ];
then
    wget https://gitlab.com/sysimm/mafft/-/archive/v7.520/mafft-v7.520.tar.gz
    tar xf mafft-v7.520.tar.gz
fi
cd src/
cp  ../mafft-v7.520/core/* .
cp Makefile Makefile.bk
sed -i -e 's/^PROGS /PROGS = dvtditr tbfast #/' -e 's/chmod /#chomd /' -e 's/sed /#sed /' -e 's/cp /#cp /' -e 's/-lm  -lpthread/-lm  ##-lpthread/' -e 's/$@ $(OBJTBFAST)/..\/build\/$@.js $(OBJTBFAST)/' -e 's/$@ $(OBJDVTDITR)/..\/build\/$@.js $(OBJDVTDITR)/'  Makefile
emmake make CC="emcc -s TOTAL_MEMORY=360MB  -s ASSERTIONS=1 $EM_FLAGS "

