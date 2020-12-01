#!/bin/bash

cd src/

# Remove large files (we'll pre-load the rest of the files as examples)
rm -rf ./test/intersect/sortAndNaming/bigTests

# Generate obj/*.o files
make clean
emmake make \
    BIN_DIR="../build/" \
    BT_LDFLAGS="--preload-file test@/bedtools2/test $(echo $EM_FLAGS)"

