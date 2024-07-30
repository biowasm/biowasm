#!/bin/bash

# Define the tarball URL
TARBALL=https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz

(
# Go to the directory containing this script
cd "$(dirname "${BASH_SOURCE[0]}")"

# Check if mummer.tar.gz exists
if [ ! -f mummer.tar.gz ]; then
    # Download the tarball to the specified directory
    wget -O mummer.tar.gz $TARBALL
    echo "Download complete."
else
    echo "mummer.tar.gz already exists. Skipping download."
fi

# Extract the tarball to src
tar -xzvf mummer.tar.gz -C src --strip-components=1
echo "Extraction complete."

# Clean up the tarball file
# rm mummer.tar.gz
)

# Create the build directory if it doesn't exist
mkdir -p ../build

# EM_FLAGS="-s USE_ZLIB=1 -s INVOKE_RUN=0 -s FORCE_FILESYSTEM=1 -s EXPORTED_RUNTIME_METHODS=["callMain","FS","PROXYFS","WORKERFS"] -s MODULARIZE=1 -s ENVIRONMENT="web,worker" -s ALLOW_MEMORY_GROWTH=1 -lworkerfs.js -lproxyfs.js"
emcc src/umd/nucmer_main.cc src/umd/nucmer.cc src/essaMEM/sparseSA.cpp src/tigr/postnuc.cc src/tigr/mgaps.cc src/essaMEM/sssort_compact.cc src/tigr/tigrinc.cc src/tigr/sw_align.cc \
    -I. -Iinclude -Iinclude/mummer \
    -o ../build/mummer.js \
    -O2 \
    $EM_FLAGS
echo "Compilation complete."
