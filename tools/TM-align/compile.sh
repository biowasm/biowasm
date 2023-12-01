# Download test data
mkdir -p testdata
cd testdata
wget -c "https://zhanggroup.org/TM-align/example/strA.pdb" "https://zhanggroup.org/TM-align/example/strB.pdb"
cd -

# Download code and compile it
wget -c "https://zhanggroup.org/TM-align/TMalign.cpp"
em++ -static -O3 -ffast-math -lm $EM_FLAGS \
  --preload-file testdata@/TMalign/testdata
  -o ../build/TMalign.js TMalign.cpp
