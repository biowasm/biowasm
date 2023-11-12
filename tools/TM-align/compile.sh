mkdir -p build src/testdata
cd src/testdata
wget -c "https://zhanggroup.org/TM-align/example/strA.pdb"
wget -c "https://zhanggroup.org/TM-align/example/strB.pdb"
cd ..
wget -c "https://zhanggroup.org/TM-align/TMalign.cpp"
em++ -static -O3 -ffast-math -lm $EM_FLAGS --preload-file  testdata@/TMalign/testdata -o ../build/TMalign.js TMalign.cpp
echo done
