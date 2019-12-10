DIR_TOOLS = tools

# Clean
clean:
	rm -rf $(DIR_TOOLS)/*/build/

all: bedtools2 bhtsne htslib samtools seqtk wgsim

init:
	@ \
	echo "——————————————————————————————————————————————————"; \
	echo "🧬 Updating git submodules..."; \
	echo "——————————————————————————————————————————————————"; \
	git submodule update --init --recursive; \
	git submodule status; \

bedtools2 bhtsne htslib samtools seqtk wgsim: init
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
