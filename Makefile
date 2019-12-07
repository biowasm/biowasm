DIR_TOOLS = tools


# Clean
# TODO: update dir_build
clean:
	rm -rf $(DIR_BUILD)

all: wgsim seqtk samtools bedtools2 bhtsne

init:
	@ \
	echo "——————————————————————————————————————————————————"; \
	echo "🧬 Updating git submodules..."; \
	echo "——————————————————————————————————————————————————"; \
	git submodule update --init --recursive; \
	git submodule status; \

wgsim seqtk: init
	@ \
	. ./shared.sh; \
	cd $(DIR_TOOLS)/$@/; \
	mkdir -p build; \
	\
	echo "\n——————————————————————————————————————————————————"; \
	echo "🧬 Applying patches..."; \
	echo "——————————————————————————————————————————————————"; \
	test -f patch && (cd src && git stash && git apply -v ../patch && cd ..) || echo "No patches"; \
	\
	echo "\n——————————————————————————————————————————————————"; \
	echo "🧬 Compiling to WebAssembly..."; \
	echo "——————————————————————————————————————————————————"; \
	./compile.sh


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
	# Patch: Reset "opt" variables so that it works properly when call main() multiple times
	# Also, use autoheader/autoconf to generate config.h.in and configure
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
	emcc $(DIR_TOOLS)/$@/samtools.o \
	  -o $(DIR_BUILD)/$@/$@.html \
	  $(EM_FLAGS) \
	  -s ERROR_ON_UNDEFINED_SYMBOLS=0


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
	emcc $(DIR_TOOLS)/$@/obj/*.o \
	  -o $(DIR_BUILD)/$@/$@.html \
	  $(EM_FLAGS) \
	  -s ERROR_ON_UNDEFINED_SYMBOLS=0


# ==============================================================================
# Tertiary analysis
# ==============================================================================

# ------------------------------------------------------------------------------
# bhtsne: t-SNE
# ------------------------------------------------------------------------------

bhtsne: init
	# Compile to .wasm and pre-load sample data
	cd $(DIR_TOOLS)/$@; \
	  sed -i 's/t-sne:/t-sne.html:/g' Makefile; \
	  gunzip brain8.snd.gz; \
	  emmake make \
		PROG="t-sne.html" \
		CC=emcc CXX=em++ \
		CFLAGS+="-s USE_ZLIB=1" \
		LIBS="-s USE_ZLIB=1 -lm --preload-file brain8.snd";

	# Move files to build folder
	for ext in data html js wasm; do \
	  mv $(DIR_TOOLS)/$@/t-sne.$$ext $(DIR_BUILD)/$@/; \
	done
