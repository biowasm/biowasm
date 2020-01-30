#!/bin/bash

cd src/

# Build:
#   - NO_TBB=1: disable Threading Building Blocks library
#   - POPCNT_CAPABILITY=0: popcnt (used in bt2_idx.h) is an assembly instruction; not supported by WebAssembly
emmake make bowtie2-align-s \
    NO_TBB=1 \
    POPCNT_CAPABILITY=0 
