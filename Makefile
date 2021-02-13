DIR_TOOLS = tools
TOOLS = $(notdir $(wildcard tools/*))
VERSION := $(VERSION)
BRANCH := $(if $(BRANCH),$(BRANCH),v$(VERSION))
TARGET := $(if $(TARGET),$(TARGET),default)

# Clean
clean:
	rm -rf $(DIR_TOOLS)/*/build/

all: ${TOOLS}

init:
	@git submodule update --init --recursive
	@git submodule status

${TOOLS}:
	@test -n "$(VERSION)"
	@./compile.sh $@ $(VERSION) $(BRANCH) $(TARGET)
