DIR_BUILD = build
DIR_TOOLS = tools


# ==============================================================================
# Initialization
# ==============================================================================

# Clean
clean:
	rm -rf $(DIR_BUILD)

# Initialize repo + folders
init:
	git submodule update --init --recursive
	for tool in $(DIR_TOOLS)/*; do \
	  mkdir -p $(DIR_BUILD)/`basename $$tool`; \
	done


# ==============================================================================
# Compile tools to WebAssembly
# ==============================================================================

# ------------------------------------------------------------------------------
# wgsim: Simulate sequencing reads
# ------------------------------------------------------------------------------

wgsim: init
	emcc $(DIR_TOOLS)/$@/$@.c \
	  -s USE_ZLIB=1 -lm \
	  -O2 -Wall \
	  -o $(DIR_BUILD)/$@/$@.js


# ------------------------------------------------------------------------------
# seqtk: FASTA/FASTQ wrangling and QC (C)
# ------------------------------------------------------------------------------

seqtk: init
	emcc $(DIR_TOOLS)/$@/$@.c \
	  -s USE_ZLIB=1 \
	  -s FORCE_FILESYSTEM=1 \
	  -s ALLOW_MEMORY_GROWTH=1 \
	  -o $(DIR_BUILD)/$@/$@.js


# ------------------------------------------------------------------------------
# samtools: SAM/BAM wrangling and QC (C)
# TODO: figure out issues with "--disable-bz2 --disable-lzma"
# ------------------------------------------------------------------------------

htslib: init
	# Install dependencies
	apt-get install -y zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev
	cd $(DIR_TOOLS)/$@/; \
	  autoheader; \
	  autoconf; \
	  emconfigure ./configure CFLAGS="-s USE_ZLIB=1" --disable-bz2 --disable-lzma

samtools: htslib
	# - Need to reset "opt" variables so that it works properly when call main() multiple times
	# - Use autoheader/autoconf to generate config.h.in and configure
	cd $(DIR_TOOLS)/$@/; \
	  sed -i "s/int ret = 0;/int ret = 0; optind = 1; opterr = 1; optopt = 0;/g" bamtk.c; \
	  autoheader; \
	  autoconf -Wno-syntax; \
	  emconfigure ./configure --without-curses CFLAGS="-s USE_ZLIB=1"; \
	  emmake make

	# Rename output to .o so it's recognizable by Emscripten
	cd $(DIR_TOOLS)/$@/; \
	  cp samtools samtools.o
	# Generate .wasm/.js files 
	emcc \
	  -o $(DIR_BUILD)/$@/$@.html $(DIR_TOOLS)/$@/samtools.o \
	  -s USE_ZLIB=1 \
	  -s ERROR_ON_UNDEFINED_SYMBOLS=0 \
	  -s INVOKE_RUN=0 \
	  -s ALLOW_MEMORY_GROWTH=1


# ------------------------------------------------------------------------------
# bedtools2: BED wrangling (C)
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
# ------------------------------------------------------------------------------

bedtools2: init
	# Build all tools
	cd $(DIR_TOOLS)/$@/; \
	  sed -i 's/^CXX.*$/CXX = emcc -s USE_ZLIB=1/' Makefile; \
	  emmake make;

	# Generate .wasm/.js files
	emcc \
	  -s USE_ZLIB=1 \
	  -s FORCE_FILESYSTEM=1 \
	  $(DIR_TOOLS)/$@/obj/*.o \
	  -o $(DIR_BUILD)/$@/$@.html \
	  -s ERROR_ON_UNDEFINED_SYMBOLS=0
