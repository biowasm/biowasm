#!/bin/bash

# ------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------

AIOLI=("1.0.0" "1.1.0" "1.1.3" "1.2.1" "1.3.0")

# Format: toolName toolVersion toolBranch
TOOLS=$(cat <<EOF
	seqtk        1.2           v1.2
	seqtk        1.3           v1.3
EOF
)

	# bedtools2    2.29.2        v2.29.2

	# bhtsne       2016.08.22    1a62a5d

	# fastp        0.20.1        0.20.1

	# samtools     1.10          1.10

	# seq-align    2017.10.18    dc41988

	# seqtk        1.2           v1.2
	# seqtk        1.3           v1.3

	# wgsim        2011.10.17    a12da33

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

	# Go to branch/tag of interest (cleanup from previous iteration)
	cd tools/${toolName}/src
	git reset --hard
	git clean -f -d
	git fetch --all
	git checkout "$toolBranch"

	# Build it
	cd ../../..
	make "$toolName"
	ls -lah tools/${toolName}/build/
	mkdir -p public/${toolName}/${toolVersion}/ public/${toolName}/latest/
	cp tools/${toolName}/build/* public/${toolName}/${toolVersion}/
	cp tools/${toolName}/build/* public/${toolName}/latest/
done <<< "$TOOLS"

# ------------------------------------------------------------------------------
# Generate CDN files for Aioli
# ------------------------------------------------------------------------------
git clone "https://github.com/biowasm/aioli.git"
cd aioli/
for version in ${AIOLI[@]};
do
	git checkout "v$version"
	dir_out="../public/aioli/$version"
	mkdir -p "$dir_out/"
	cp aioli{,.worker}.js "$dir_out/"
done
dir_out="../public/aioli/latest"
mkdir -p "$dir_out/"
cp aioli{,.worker}.js "$dir_out/"
cd ../


# ------------------------------------------------------------------------------
# Generate index
# ------------------------------------------------------------------------------
cd public/
( tree --du -h | grep -v -E "index.html|404.html|.ico" | tail +2 ) > index.txt
