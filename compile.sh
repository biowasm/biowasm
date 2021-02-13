#!/bin/bash

# This script compiles a tool to WebAssembly and is invoked by the Makefile.

TOOL=${1?Usage: ./compile.sh toolName targetName}
TARGET=${2?Usage: ./compile.sh toolName targetName}
DIR_TOOLS=tools/

# Load target Emscripten flags
. ./config/shared.$TARGET.sh

# Apply patches, if any
cd "${DIR_TOOLS}/${TOOL}/"
mkdir -p build
echo "——————————————————————————————————————————————————"
echo "🧬 Applying patches..."
echo "——————————————————————————————————————————————————"
test -f patch && (cd src && git stash && git apply -v ../patch && cd ..) || echo "No patches"

# Launch tool's compile.sh script
echo "——————————————————————————————————————————————————"
echo "🧬 Compiling to WebAssembly..."
echo "——————————————————————————————————————————————————"
echo $EM_FLAGS
./compile.sh

# Generate glue code for all .js files. The default thisProgram is "./this.program"
# so we update it to match the tool name for convenience
SHARED_JS_EXTRA=$(cat <<EOF
    if(typeof Module["thisProgram"] == "undefined")
        Module["thisProgram"] = "./${TOOL}";
EOF
)

for glueCode in build/*.js;
do
    echo "$SHARED_JS_EXTRA" > $glueCode.tmp
    cat ../../config/shared.js $glueCode.tmp $glueCode > $glueCode.final
    mv $glueCode.final $glueCode
    rm $glueCode.tmp
done
