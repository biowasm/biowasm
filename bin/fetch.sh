#!/bin/bash
# Clone a tool's source at a given branch/tag. Clone URLs live in biowasm.json.
set -euo pipefail

usage="Usage: bin/fetch.sh <dir> <url> <ref>"
DIR=${1?$usage}
URL=${2?$usage}
REF=${3?$usage}

# Some GNU repos have submodules that point to https instead of git, which seems to be more flaky.
# Submodules wait until after checkout, since the default branch can pin different ones than the tag.
REWRITE="-c url.git://git.sv.gnu.org/.insteadOf=https://git.sv.gnu.org/git/"

if [[ ! -e "$DIR/.git" ]]; then
	rm -rf "$DIR"
	git clone "$URL" "$DIR"
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
git -C "$DIR" $REWRITE submodule update --init --recursive --force
