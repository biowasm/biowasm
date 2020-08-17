#!/bin/bash

TOOLS=$(cat <<EOF
	seqtk        1.2           v1.2
	seqtk        1.3           v1.3

	fastp        0.20.1        0.20.1

	seq-align    2017.10.18    dc41988

	bhtsne       2016.08.22    1a62a5d

	samtools     1.10          1.10

	bedtools2    2.29.2        v2.29.2
EOF
)

# Setup repos and dependencies
make init
sudo apt-get install -y tree

# Compile each tool
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

# Generate index
cd public/
(
	tree | \
		grep -v -E "index|index.html|.ico" | \
		tail +2; \
) > index
