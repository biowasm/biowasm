#!/bin/bash

# Compile lastz to WebAssembly
# Common EM_FLAGS are provided by biowasm compile.py
# This script runs from tools/lastz/src/ directory

# Get version info
VERSION_MAJOR=$(grep "^VERSION_MAJOR" src/version.mak | cut -d= -f2)
VERSION_MINOR=$(grep "^VERSION_MINOR" src/version.mak | cut -d= -f2)
VERSION_SUBMINOR=$(grep "^VERSION_SUBMINOR" src/version.mak | cut -d= -f2)
REVISION_DATE=$(grep "^REVISION_DATE" src/version.mak | cut -d= -f2)

echo "Compiling lastz version $VERSION_MAJOR.$VERSION_MINOR.$VERSION_SUBMINOR"

# Version flags
VERSION_FLAGS="-DVERSION_MAJOR=\"$VERSION_MAJOR\" -DVERSION_MINOR=\"$VERSION_MINOR\" -DVERSION_SUBMINOR=\"$VERSION_SUBMINOR\" -DREVISION_DATE=\"$REVISION_DATE\" -DSUBVERSION_REV=\"\""

# Change to source directory
cd src

# Compile all source files
emcc -O3 -c \
  -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE \
  -Dscore_type=I \
  $VERSION_FLAGS \
  lastz.c infer_scores.c seeds.c pos_table.c quantum.c seed_search.c diag_hash.c \
  chain.c gapped_extend.c tweener.c masking.c segment.c edit_script.c \
  identity_dist.c coverage_dist.c continuity_dist.c \
  output.c gfa.c lav.c axt.c maf.c cigar.c sam.c genpaf.c text_align.c align_diffs.c \
  utilities.c dna_utilities.c sequences.c capsule.c

# Create build directory and link
mkdir -p ../../build
emcc -O3 *.o -o ../../build/lastz.js $EM_FLAGS

echo "lastz compiled successfully!"




