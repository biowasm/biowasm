#!/bin/bash

# This script compiles the bioinformatics tools to WebAssembly and is invoked by
# the GitHub Actions pipelines in `.github/workflows/`.

# ------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------

DIR_CDN="cloudflare/cdn/public"

AIOLI=("1.4.1")

# Format: toolName toolVersion toolBranch
TOOLS=$(cat <<EOF
	bedtools2    2.29.2        v2.29.2

	bhtsne       2016.08.22    1a62a5d

	bowtie2      2.4.2         v2.4.2

	fastp        0.20.1        0.20.1

	samtools     1.10          1.10

	seq-align    2017.10.18    dc41988

	seqtk        1.2           v1.2
	seqtk        1.3           v1.3

	ssw          1.2.4         ad452ea

	wgsim        2011.10.17    a12da33
EOF
)

# ------------------------------------------------------------------------------
# Setup repos and dependencies
# ------------------------------------------------------------------------------
make init
sudo apt-get install -y tree

# ------------------------------------------------------------------------------
# Compile each tool
# ------------------------------------------------------------------------------
while read toolName toolVersion toolBranch;
do
	if [[ "$toolName" == "" ]]; then
		continue;
	fi

	echo "================================================================"
	echo "Processing $toolName/$toolVersion @ $toolBranch"
	echo "================================================================"

	# Build it
	VERSION="$toolVersion" BRANCH="$toolBranch" make "$toolName"

	# Copy over config.json file if it exists
	toolConfig="tools/${toolName}/configs/${toolVersion}.json"
	[[ -f "$toolConfig" ]] && cp "$toolConfig" tools/${toolName}/build/config.json

	# Copy files over to the expected CDN folder
	ls -lah tools/${toolName}/build/
	mkdir -p ${DIR_CDN}/${toolName}/${toolVersion}/ ${DIR_CDN}/${toolName}/latest/
	cp tools/${toolName}/build/* ${DIR_CDN}/${toolName}/${toolVersion}/
	cp tools/${toolName}/build/* ${DIR_CDN}/${toolName}/latest/
done <<< "$TOOLS"

# ------------------------------------------------------------------------------
# Generate CDN files for Aioli
# ------------------------------------------------------------------------------
git clone "https://github.com/biowasm/aioli.git"
cd aioli/
for version in ${AIOLI[@]};
do
	git checkout "v$version"
	dir_out="../$DIR_CDN/aioli/$version"
	mkdir -p "$dir_out/"
	cp aioli{,.worker}.js "$dir_out/"
done
dir_out="../$DIR_CDN/aioli/latest"
mkdir -p "$dir_out/"
cp aioli{,.worker}.js "$dir_out/"
cd ../

# ------------------------------------------------------------------------------
# Generate index
# ------------------------------------------------------------------------------
cd "$DIR_CDN"
( echo "cdn.biowasm.com"; date; tree --charset=ascii --du -h | grep -v -E "index.html|index|404.html|.ico" | tail +2 ) > index
