#!/bin/bash

# TODO: look into:
#   shared:WARNING: emcc: cannot find library "bz2"
#   shared:WARNING: emcc: cannot find library "lzma"
#   warning: undefined symbol: bam_aux_append
#   warning: undefined symbol: bam_aux_get
#   warning: undefined symbol: bam_copy1
#   warning: undefined symbol: bam_endpos
#   warning: undefined symbol: bam_hdr_destroy
#   warning: undefined symbol: bam_hdr_init
#   warning: undefined symbol: bgzf_hopen
#   warning: undefined symbol: bgzf_read
#   warning: undefined symbol: cram_get_refs
#   warning: undefined symbol: hopen_callback
#   warning: undefined symbol: hts_close
#   warning: undefined symbol: hts_idx_destroy
#   warning: undefined symbol: hts_itr_destroy
#   warning: undefined symbol: hts_itr_next
#   warning: undefined symbol: hts_open
#   warning: undefined symbol: hts_open_callback
#   warning: undefined symbol: hts_set_fai_filename
#   warning: undefined symbol: hts_set_opt
#   warning: undefined symbol: sam_hdr_read
#   warning: undefined symbol: sam_hdr_write
#   warning: undefined symbol: sam_index_load
#   warning: undefined symbol: sam_itr_queryi
#   warning: undefined symbol: sam_read1
#   warning: undefined symbol: sam_write1

cd src/

# Patch: Reset "opt" variables so that it works properly when call main() multiple times
# Also, use autoheader/autoconf to generate config.h.in and configure
autoheader
autoconf -Wno-syntax
emconfigure ./configure --without-curses --with-htslib="../../htslib/src/" CFLAGS="-s USE_ZLIB=1"
emmake make CC=emcc AR=emar

# Rename output to .o so it's recognizable by Emscripten
cp samtools samtools.o

# Generate .wasm/.js files 
emcc samtools.o \
    -o ../build/samtools.html \
    $EM_FLAGS \
    -s ERROR_ON_UNDEFINED_SYMBOLS=0
