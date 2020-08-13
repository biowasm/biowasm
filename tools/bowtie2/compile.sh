#!/bin/bash

cd src/

# Cloudflare currently has a 10MB/file limit so remove unneeded files
rm example/reads/combined_reads.bam
rm example/reads/longreads.fq

cat example/reads/reads_1.fq | head -n4000 | gzip > example/reads/reads_1.fq.tmp
cat example/reads/reads_2.fq | head -n4000 | gzip > example/reads/reads_2.fq.tmp
mv example/reads/reads_1.fq.tmp example/reads/reads_1.fq
mv example/reads/reads_2.fq.tmp example/reads/reads_2.fq

# Build:
#   - NO_TBB=1: disable Threading Building Blocks library
#   - POPCNT_CAPABILITY=0: popcnt (used in bt2_idx.h) is an assembly instruction; not supported by WebAssembly
emmake make bowtie2-align-s \
    NO_TBB=1 \
    POPCNT_CAPABILITY=0 
