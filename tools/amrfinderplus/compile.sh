#!/bin/bash

set -e

# amrfinder and stxtyper drive these through system(), which WebAssembly has no equivalent of,
# so the programs they drive are built instead and the pipeline is run from JavaScript
PROGRAMS="amr_report dna_mutation fasta_check fasta_extract fasta2parts gff_check"

# The Makefile passes the version in as a macro; it is not read from version.txt
VERSION=$(cat version.txt)

for program in $PROGRAMS; do
	objs=$(grep "^${program}OBJS=" Makefile | cut -d= -f2)
	em++ -O2 -std=gnu++20 -DNDEBUG -fexceptions -DSVN_REV=\"$VERSION\" \
		${objs//.o/.cpp} \
		-o ../build/$program.js \
		$EM_FLAGS -sERROR_ON_UNDEFINED_SYMBOLS=0
done
