
## below is only for muscle v5.1.0
#cd src/
#git apply ../patch/5.1.0.patch
#emmake make  CC=emcc CXX=em++


## below is only for muscle v3.8.1551
# Defines the env variable $EM_FLAGS
source ../../bin/shared.sh
mkdir -p src build
wget https://drive5.com/muscle/muscle_src_3.8.1551.tar.gz
tar xf muscle_src_3.8.1551.tar.gz -C src
sed -i -r -e 's/strip muscle/#strip muscle/' -e 's/ -o muscle / -o ..\/build\/muscle.js \$\(EM_FLAGS\) /' src/Makefile
cd src
emmake make GPP=emcc EM_FLAGS="$EM_FLAGS"
cd ..
rm muscle_src_3.8.1551.tar.gz

