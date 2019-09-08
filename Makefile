
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
# TODO: figure out issues with "--disable-bz2 --disable-lzma"
# ------------------------------------------------------------------------------

htslib: init
	@# Install dependencies
	apt-get install -y zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev
	cd $(DIR_TOOLS)/$@/; \
	  autoheader; \
	  autoconf; \
	  emconfigure ./configure CFLAGS="-s USE_ZLIB=1" --disable-bz2 --disable-lzma

samtools: htslib
	@# - Need to reset "opt" variables so that it works properly when call main() multiple times
	@# - Use autoheader/autoconf to generate config.h.in and configure
	cd $(DIR_TOOLS)/$@/; \
	  sed -i "s/int ret = 0;/int ret = 0; optind = 1; opterr = 1; optopt = 0;/g" bamtk.c; \
	  autoheader; \
	  autoconf -Wno-syntax; \
	  emconfigure ./configure --without-curses CFLAGS="-s USE_ZLIB=1"; \
	  emmake make

	@# Rename output to .o so it's recognizable by Emscripten
	cd $(DIR_TOOLS)/$@/; \
	  cp samtools samtools.o
	@# Generate .wasm/.js files 
	emcc \
	  -o $(DIR_BUILD)/$@/$@.html $(DIR_TOOLS)/$@/samtools.o \
	  -s USE_ZLIB=1 \
	  -s ERROR_ON_UNDEFINED_SYMBOLS=0 \
	  -s INVOKE_RUN=0 \
	  -s ALLOW_MEMORY_GROWTH=1
