#!/bin/bash
# Clone a tool's source at a given branch/tag. Clone URLs live in biowasm.json.
set -euo pipefail

usage="Usage: bin/fetch.sh <dir> <url> <ref>"
DIR=${1?$usage}
URL=${2?$usage}
REF=${3?$usage}

# Blobless clones are much smaller, but not every host supports them
if [[ ! -e "$DIR/.git" ]]; then
	rm -rf "$DIR"
	git clone --filter=blob:none --recurse-submodules "$URL" "$DIR" ||
		git clone --recurse-submodules "$URL" "$DIR"
fi

# Without its own .git, the commands below would act on the biowasm repo instead
if [[ "$(git -C "$DIR" rev-parse --show-toplevel)" != "$(cd "$DIR" && pwd)" ]]; then
	echo "ERROR: '$DIR' is not its own git repo; refusing to overwrite it." >&2
	exit 1
fi

# Move to the requested ref and drop anything the last build left behind
git -C "$DIR" fetch --tags --force origin
git -C "$DIR" checkout --force "$REF"
git -C "$DIR" clean -xdff
git -C "$DIR" submodule update --init --recursive
