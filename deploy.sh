#!/bin/bash

# This script compiles the bioinformatics tools to WebAssembly and is invoked by
# the GitHub Actions pipelines in `.github/workflows/`.

# ------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------

DIR_TOOLS="config/tools.json"
DIR_CDN="cloudflare/cdn/public"
URL_CDN="https://cdn.biowasm.com/v2"

# ------------------------------------------------------------------------------
# Setup repos and dependencies
# ------------------------------------------------------------------------------
[[ "$CACHE_DISABLED" == "true" ]] && make init
sudo apt-get install -y tree jq

# ------------------------------------------------------------------------------
# Compile each tool
# ------------------------------------------------------------------------------
echo "Running with ENV=${ENV}..."
echo "Running with CACHE_DISABLED=${CACHE_DISABLED}..."

# Load info about each tool into an array
allTools=($(jq -rc '.tools[]' $DIR_TOOLS))

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
	if [[ "$CACHE_DISABLED" == "true" ]]; then
		VERSION="$toolVersion" BRANCH="$toolBranch" make "$toolName"
	else
		[[ "$ENV" == "prd" ]] && url=$URL_CDN || url="${URL_CDN//cdn/cdn-stg}"
		for program in "${toolPrograms[@]}"; do
			curl -o tools/${toolName}/build/config.json "${url}/${toolName}/${toolVersion}/config.json"
			curl -o tools/${toolName}/build/${program}.js "${url}/${toolName}/${toolVersion}/${program}.js"
			curl -o tools/${toolName}/build/${program}.wasm "${url}/${toolName}/${toolVersion}/${program}.wasm"
			curl --fail -o tools/${toolName}/build/${program}.data "${url}/${toolName}/${toolVersion}/${program}.data"  # ignore .data failures since not all tools have .data files
		done
	fi

	# Copy files over to the expected CDN folder
	ls -lah tools/${toolName}/build/
	mkdir -p ${DIR_CDN}/${toolName}/${toolVersion}/
	cp tools/${toolName}/build/* ${DIR_CDN}/${toolName}/${toolVersion}/
done

# ------------------------------------------------------------------------------
# Generate CDN files for Aioli
# ------------------------------------------------------------------------------
git clone "https://github.com/biowasm/aioli.git"
cd aioli/

allAiolis=($(jq -rc '.aioli[]' $DIR_TOOLS))
for aioli in ${allAiolis[@]};
do
	aioliVersion=$(jq -rc '.version' <<< $tool)
	aioliBranch=$(jq -rc '.branch' <<< $tool)

	git checkout "$aioliBranch"
	dir_out="../$DIR_CDN/aioli/$aioliVersion"
	mkdir -p "$dir_out/"

	npm install
	npm run build

	cp dist/aioli{,.worker}.js "$dir_out/"
done
dir_out="../$DIR_CDN/aioli/latest"
mkdir -p "$dir_out/"
cp aioli{,.worker}.js "$dir_out/"
cd ../

# ------------------------------------------------------------------------------
# Generate index
# ------------------------------------------------------------------------------
cd "$DIR_CDN"
( echo "cdn.biowasm.com/v2/"; date; tree --charset=ascii --du -h | grep -v -E "index.html|index|404.html|.ico" | tail +2 ) > index
