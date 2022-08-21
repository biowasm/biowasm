#!/bin/bash

cd src/

# Make example tree
echo -e "((A:0.12,D:0.12):0.3,(B:0.3,C:0.5):0.4);" > tree.newick
echo -e "5\t\t\t# number of temporal constraints\nA 1999.2\t\t# the date of A is 1999.2\nB 2000.1\t\t# the date of B is 2000.1\nC l(1990.5)\t\t# the date of C is >= 1990.5 (more recent than 1990.5)\nD b(1998.21,2000.5)\t# the date of D is between 1998.21 and 2000.5\nmrca(A,B,C) u(1980)\t# the date of the most recent ancestor of A,B, and C is <= 1980 (older than 1980)" > tree.date

# All .cpp files to build
TO_BUILD=$(ls *.cpp | sed "s/.cpp/.o/")

# Build .cpp files
make clean  # repo contains precompiled .o
emmake make $TO_BUILD CXX=em++ CXXFLAGS="-O2"

# Build .js/.wasm
em++ -o ../../build/lsd2.js $TO_BUILD $EM_FLAGS -O2 \
	--preload-file tree.newick@/lsd2/tree.newick \
	--preload-file tree.date@/lsd2/tree.date
