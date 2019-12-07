DIR_TOOLS = tools

# Clean
clean:
	rm -rf tools/*/build/

all: wgsim seqtk samtools bedtools2 bhtsne

init:
	@ \
	echo "——————————————————————————————————————————————————"; \
	echo "🧬 Updating git submodules..."; \
	echo "——————————————————————————————————————————————————"; \
	git submodule update --init --recursive; \
	git submodule status; \

wgsim seqtk bhtsne bedtools2 htslib samtools: init
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
	./compile.sh; \
	\
	cd src; \
	rm a.out{,.js,.wasm}; \
	git stash;
