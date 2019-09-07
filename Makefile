
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
# seqtk: FASTA/FASTQ wrangling and QC (C)
# ------------------------------------------------------------------------------

seqtk: init
	emcc $(DIR_TOOLS)/$@/$@.c \
	  -s USE_ZLIB=1 \
	  -s FORCE_FILESYSTEM=1 \
	  -s EXTRA_EXPORTED_RUNTIME_METHODS='["ccall", "cwrap"]' \
	  -s ALLOW_MEMORY_GROWTH=1 \
	  -o $(DIR_BUILD)/$@/$@.js


# ------------------------------------------------------------------------------
# samtools: SAM/BAM wrangling and QC (C)
# ------------------------------------------------------------------------------

htslib: init
	

samtools: init
	# Dependencies
	apt-get install -y zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev

	# Build htslib
	cd $(DIR_TOOLS)/$@/

	cd htslib/
	autoheader
	autoconf
	emconfigure ./configure CFLAGS="-s USE_ZLIB=1" --disable-bz2 --disable-lzma
	# TODO: is this step needed?
	emmake make


