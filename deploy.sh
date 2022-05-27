#!/bin/bash

# This script compiles the bioinformatics tools to WebAssembly and is invoked by
# the GitHub Actions pipelines in `.github/workflows/`.

# ------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------

DIR_TOOLS="config/tools.json"
DIR_CDN="cloudflare/cdn/public"
URL_CDN="https://cdn.biowasm.com/v2"
# ENV               # Either "stg" or "prd"
# TOOLS_TO_COMPILE  # Comma-separated list of tools to recompile; the rest will use CDN as cache ('all', 'none')

# ------------------------------------------------------------------------------
# Setup repos and dependencies
# ------------------------------------------------------------------------------
[[ "$TOOLS_TO_COMPILE" == "none" ]] && TOOLS_TO_COMPILE=""
sudo apt-get install -y tree jq

# ------------------------------------------------------------------------------
# Compile each tool
# ------------------------------------------------------------------------------
echo "Running with ENV=${ENV}..."
echo "Running with TOOLS_TO_COMPILE=${TOOLS_TO_COMPILE}..."

# Load info about each tool into an array
allTools=($(jq -rc '.tools[]' $DIR_TOOLS))
# Load list of tools to compile into an array
IFS="," read -r -a TOOLS_TO_COMPILE <<< "$TOOLS_TO_COMPILE"

# Initialize repo for tools of interest
for tool in "${TOOLS_TO_COMPILE[@]}"; do
	make init $tool
done

# Build each tool
for tool in "${allTools[@]}";
do
	# Parse tool info
	toolName=$(jq -rc '.name' <<< $tool)
	toolVersion=$(jq -rc '.version' <<< $tool)
	toolBranch=$(jq -rc '.branch' <<< $tool)
	toolPrograms=$(jq -rc '.programs' <<< $tool)
	[[ "$toolPrograms" == "null" ]] && toolPrograms="[\"$toolName\"]"
	toolPrograms=($(jq -rc '.[]' <<< $toolPrograms))

	# Compile it to WebAssembly or fetch pre-compiled from existing CDN!
	if ( [[ "${TOOLS_TO_COMPILE[0]}" == "all" ]] || [[ " ${TOOLS_TO_COMPILE[@]} " =~ " ${toolName} " ]] ); then
		VERSION="$toolVersion" BRANCH="$toolBranch" make "$toolName"
	else
		mkdir -p tools/${toolName}/build/
		[[ "$ENV" == "prd" ]] && url=$URL_CDN || url="${URL_CDN//cdn/cdn-stg}"
		curl -s -o tools/${toolName}/build/config.json "${url}/${toolName}/${toolVersion}/config.json"
		for program in "${toolPrograms[@]}"; do
			curl -s -o tools/${toolName}/build/${program}.js "${url}/${toolName}/${toolVersion}/${program}.js"
			curl -s -o tools/${toolName}/build/${program}.wasm "${url}/${toolName}/${toolVersion}/${program}.wasm"
			curl -s --fail -o tools/${toolName}/build/${program}.data "${url}/${toolName}/${toolVersion}/${program}.data"  # ignore .data failures since not all tools have .data files
		done
	fi
	echo "> tools/${toolName}/build/"
	ls -lah tools/${toolName}/build/

	# Copy files over to the expected CDN folder
	mkdir -p ${DIR_CDN}/${toolName}/${toolVersion}/
	cp tools/${toolName}/build/* ${DIR_CDN}/${toolName}/${toolVersion}/
	echo "> ${DIR_CDN}/${toolName}/${toolVersion}/"
	ls -lah ${DIR_CDN}/${toolName}/${toolVersion}/
done

# ------------------------------------------------------------------------------
# Generate CDN files for Aioli
# ------------------------------------------------------------------------------
allAiolis=($(jq -rc '.aioli[]' $DIR_TOOLS))

git clone "https://github.com/biowasm/aioli.git"
cd aioli/
dir_out_latest="../$DIR_CDN/aioli/latest"
mkdir -p "$dir_out_latest/"

for aioli in ${allAiolis[@]};
do
	aioliVersion=$(jq -rc '.version' <<< $aioli)
	aioliBranch=$(jq -rc '.branch' <<< $aioli)
	aioliLatest=$(jq -rc '.latest' <<< $aioli)

	git checkout "$aioliBranch"
	dir_out="../$DIR_CDN/aioli/$aioliVersion"
	mkdir -p "$dir_out/"

	npm install
	npm run build

	cp dist/aioli{,.worker}.js "$dir_out/"
	if [[ "$aioliLatest" != "false" ]]; then
		cp dist/aioli{,.worker}.js "$dir_out_latest/"
	fi
done
cd ../


# ------------------------------------------------------------------------------
# Generate index
# ------------------------------------------------------------------------------
cd "$DIR_CDN"
prefix="cdn"; [[ "$ENV" == "stg" ]] && prefix="cdn-stg"
( echo "$prefix.biowasm.com/v2/"; date; tree --charset=ascii --du -h | grep -v -E "index.html|index|404.html|.ico" | tail +2 ) > index
