#!/bin/bash

cd src/

# Remove unneeded files to keep .data files small so they load quickly
rm example/reads/combined_reads.bam
rm example/reads/longreads.fq
rm -rf example/index

cat example/reads/reads_1.fq | head -n100 > example/reads/reads_1.fq.tmp
cat example/reads/reads_2.fq | head -n100 > example/reads/reads_2.fq.tmp
mv example/reads/reads_1.fq.tmp example/reads/reads_1.fq
mv example/reads/reads_2.fq.tmp example/reads/reads_2.fq

# Build:
#   - NO_TBB=1: disable Threading Building Blocks library
#   - POPCNT_CAPABILITY=0: popcnt (used in bt2_idx.h) is an assembly instruction; not supported by WebAssembly
emmake make bowtie2-align-s \
    NO_TBB=1 \
    POPCNT_CAPABILITY=0 \
    WASM_FLAGS="$(echo $EM_FLAGS) --preload-file example/@/bowtie2/example/ -s ASSERTIONS=1 -s DISABLE_EXCEPTION_CATCHING=0"
