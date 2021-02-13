#!/bin/bash

# This script compiles a tool to WebAssembly and is invoked by the Makefile.

usage="Usage: ./compile.sh toolName toolVersion toolBranch targetName"
TOOL=${1?$usage}
VERSION=${2?$usage}
BRANCH=${3?$usage}
TARGET=${4?$usage}
DIR_TOOLS=tools/

# ------------------------------------------------------------------------------
# Initialize
# ------------------------------------------------------------------------------

# Load target Emscripten flags
. ./config/shared.$TARGET.sh

# Prep build/ folder
cd "${DIR_TOOLS}/${TOOL}/"
mkdir -p build

# ------------------------------------------------------------------------------
# Setup codebase
# ------------------------------------------------------------------------------

echo "——————————————————————————————————————————————————"
echo "🧬 Processing $TOOL v$VERSION @ branch '$BRANCH'..."
echo "——————————————————————————————————————————————————"

cd src/

# Go to branch/tag of interest (clean up previous iterations)
git reset --hard
git clean -f -d
git fetch --all
git checkout "$BRANCH"

# Apply patches, if any
patch_file=../patches/$VERSION
if [[ -f "$patch_file" ]]; then
    echo "Applying patch file <$patch_file>"
    git apply -v $patch_file
else
    echo "No patch file found"
fi

cd ../

# ------------------------------------------------------------------------------
# Compile tool
# ------------------------------------------------------------------------------

echo "——————————————————————————————————————————————————"
echo "🧬 Compiling $TOOL v$VERSION to WebAssembly..."
echo "——————————————————————————————————————————————————"
./compile.sh

# ------------------------------------------------------------------------------
# Generate glue code for all .js files
# ------------------------------------------------------------------------------

# Default thisProgram == "./this.program"; update it to match tool name for convenience
cat <<EOF > glueCode.extra
    if(typeof Module["thisProgram"] == "undefined")
        Module["thisProgram"] = "./${TOOL}";
EOF

# Combine all glue code into 1 .js file
for glueCode in build/*.js; do
    cat ../../config/shared.js glueCode.extra $glueCode > $glueCode.final
    mv $glueCode.final $glueCode
done
rm glueCode.extra
