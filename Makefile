DIR_TOOLS = tools
TOOLS = $(notdir $(wildcard tools/*))
VERSION := $(VERSION)
BRANCH := $(if $(BRANCH),$(BRANCH),v$(VERSION))
TARGET := $(if $(TARGET),$(TARGET),default)
TOOL := $(if $(TOOL),tools/$(TOOL)/src/,)

# Clean
clean:
	rm -rf $(DIR_TOOLS)/*/build/

all: ${TOOLS}

init:
	@git submodule update --init --recursive $(TOOL)
	@git submodule status $(TOOL)

${TOOLS}:
	@test -n "$(VERSION)"
	@./compile.sh $@ $(VERSION) $(BRANCH) $(TARGET)
