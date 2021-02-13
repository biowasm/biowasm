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
	./compile.sh $@ $(TARGET)
