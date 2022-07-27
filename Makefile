DIR_TOOLS = tools
TOOL := $(if $(TOOL),tools/$(TOOL)/src/,)

help:
	@echo "Run 'make init <toolname>' to initialize the repo for a tool of interest"

clean:
	rm -rf $(DIR_TOOLS)/*/build/

init:
	@git submodule update --init --recursive $(TOOL)
	@git submodule status $(TOOL)
