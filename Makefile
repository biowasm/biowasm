DIR_TOOLS = tools
TOOLS = $(notdir $(wildcard tools/*))
TARGET := $(if $(TARGET),$(TARGET),default)

# Clean
clean:
	rm -rf $(DIR_TOOLS)/*/build/

all: ${TOOLS}

init:
	@ \
	echo "——————————————————————————————————————————————————"; \
	echo "🧬 Updating git submodules..."; \
	echo "——————————————————————————————————————————————————"; \
	git submodule update --init --recursive; \
	git submodule status; \

${TOOLS}:
	@ \
	. ./config/shared.$(TARGET).sh; \
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
	for glueCode in build/*.js; do \
		cat ../../config/shared.js $$glueCode > $$glueCode.tmp; \
		mv $$glueCode.tmp $$glueCode; \
	done
