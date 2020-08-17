#!/bin/bash

TOOLS=$(cat <<EOF
	seqtk    1.2    v1.2
	seqtk    1.3    v1.3
EOF
)

# Setup repos and dependencies
make init
sudo apt-get install -y tree liblzma-dev

# Compile each tool
while read toolName toolVersion toolBranch;
do
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
