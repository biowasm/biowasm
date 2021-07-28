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
	toolPath=$(jq -rc '.path' <<< $tool)  # e.g. bedtools source is in bedtools2/
	[[ "$toolPath" == "null" ]] && toolPath="${toolName}"
	toolVersion=$(jq -rc '.version' <<< $tool)
	toolBranch=$(jq -rc '.branch' <<< $tool)
	toolPrograms=$(jq -rc '.programs' <<< $tool)
	[[ "$toolPrograms" == "null" ]] && toolPrograms="[\"$toolName\"]"
	toolPrograms=($(jq -rc '.[]' <<< $toolPrograms))

	# Compile it to WebAssembly or fetch pre-compiled from existing CDN!
	if [[ "$CACHE_DISABLED" == "true" ]]; then
		VERSION="$toolVersion" BRANCH="$toolBranch" make "$toolPath"
	else
		mkdir -p tools/${toolPath}/build/
		[[ "$ENV" == "prd" ]] && url=$URL_CDN || url="${URL_CDN//cdn/cdn-stg}"
		curl -s -o tools/${toolPath}/build/config.json "${url}/${toolName}/${toolVersion}/config.json"
		for program in "${toolPrograms[@]}"; do
			curl -s -o tools/${toolPath}/build/${program}.js "${url}/${toolName}/${toolVersion}/${program}.js"
			curl -s -o tools/${toolPath}/build/${program}.wasm "${url}/${toolName}/${toolVersion}/${program}.wasm"
			curl -s --fail -o tools/${toolPath}/build/${program}.data "${url}/${toolName}/${toolVersion}/${program}.data"  # ignore .data failures since not all tools have .data files
		done
	fi

	echo "> tools/${toolPath}/build/"
	ls -lah tools/${toolPath}/build/

	# Copy files over to the expected CDN folder
	mkdir -p ${DIR_CDN}/${toolName}/${toolVersion}/
	cp tools/${toolPath}/build/config.json ${DIR_CDN}/${toolName}/${toolVersion}/config.json
	# Some tools have multiple programs (e.g. ssw has smith_waterman, needleman_wunsch, and lcs)
	for program in "${toolPrograms[@]}"; do
		cp tools/${toolPath}/build/${program}.js ${DIR_CDN}/${toolName}/${toolVersion}/${program}.js
		cp tools/${toolPath}/build/${program}.wasm ${DIR_CDN}/${toolName}/${toolVersion}/${program}.wasm
		if [[ -f "tools/${toolPath}/build/${program}.data" ]]; then
			cp tools/${toolPath}/build/${program}.data ${DIR_CDN}/${toolName}/${toolVersion}/${program}.data
		fi
	done

	echo "> ${DIR_CDN}/${toolName}/${toolVersion}/"
	ls -lah ${DIR_CDN}/${toolName}/${toolVersion}/
done

# ------------------------------------------------------------------------------
# Generate CDN files for Aioli
# ------------------------------------------------------------------------------
allAiolis=($(jq -rc '.aioli[]' $DIR_TOOLS))

git clone "https://github.com/biowasm/aioli.git"
cd aioli/

for aioli in ${allAiolis[@]};
do
	aioliVersion=$(jq -rc '.version' <<< $aioli)
	aioliBranch=$(jq -rc '.branch' <<< $aioli)

	git checkout "$aioliBranch"
	dir_out="../$DIR_CDN/aioli/$aioliVersion"
	mkdir -p "$dir_out/"

	npm install
	npm run build

	cp dist/aioli{,.worker}.js "$dir_out/"
done
dir_out="../$DIR_CDN/aioli/latest"
mkdir -p "$dir_out/"
cp dist/aioli{,.worker}.js "$dir_out/"
cd ../


# ------------------------------------------------------------------------------
# Generate index
# ------------------------------------------------------------------------------
cd "$DIR_CDN"
( echo "cdn.biowasm.com/v2/"; date; tree --charset=ascii --du -h | grep -v -E "index.html|index|404.html|.ico" | tail +2 ) > index
