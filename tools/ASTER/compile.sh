[[ -e data ]] && rm data/*
mkdir -p data
cp -r example/multitree_genename.map example/multitree_genename.nw example/genetrees.tre_1.fas -t data/
make -f Makefile.emcc LDFLAGS="$EM_FLAGS --preload-file data@/ASTER/example"
mv bin/* ../build/
