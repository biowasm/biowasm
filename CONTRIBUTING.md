# Contributing

## Development Setup

Tools listed in biowasm were compiled to WebAssembly using `Emscripten 6.0.13`. Here is how to set up your dev environment:

```bash
# Fetch Emscripten docker image
docker build -t biowasm-dev -f Dockerfile-dev .

# Create the container and mount ~/wasm to /src in the container
docker run \
    -it -d \
    -p 80:80 \
    --name wasm \
    --volume ~/wasm:/src \
    biowasm-dev
```

This creates a simple CORS enabled development server serving the local `~/wasm` directory, 
ideally the parent directory of your cloned `biowasm` repo.

## Compile an existing biowasm tool

```bash
# Go into your container
docker exec -it wasm bash

# Compile seqtk
cd biowasm/
bin/compile.py --tools seqtk --versions 1.2

# This will create tools/<tool name>/build with .js/.wasm files
ls tools/seqtk/build
```


## Add a new tool to biowasm

Add a new tool to biowasm.json (we no longer use `.gitmodules` because it doesn't scale and we typically compile one or a few tools at a time):

```json
{
	"name": "minimap2",
	"url": "https://github.com/lh3/minimap2.git",
	"programs": ["minimap2"],
	"description": "Sequence alignment program",
	"versions": [
		{ "version": "2.22", "branch": "v2.22" }
	]
}
```

Then fetch the source and scaffold the build:

```bash
TOOL=minimap2
BRANCH=v2.22
REPO=https://github.com/lh3/minimap2.git

mkdir -p tools/$TOOL/examples

# Clone the tool at that tag (same thing bin/compile.py does before each build)
bin/fetch.sh tools/$TOOL/src $REPO $BRANCH

# Then create a compile.sh script that compiles the tool to WebAssembly
echo "# TODO" > tools/$TOOL/compile.sh
chmod +x tools/$TOOL/compile.sh

mkdir -p tools/$TOOL/examples
echo "#TODO" > tools/$TOOL/examples/$BRANCH.html
```

You should also modify/create the following files:

```bash
biowasm.json         Add the new tool here (with its clone `url`) so it gets deployed to the CDN

tools/<tool>/
    compile.sh       Script that will run to compile the tool to WebAssembly (use `$EM_FLAGS` for common flags)
    patches/    
        <tag>.patch  Patch applied to the code to compile it to WebAssembly; branch- or tag-specific (optional)
    examples/
        <tag>.html   Example HTML code for running the tool; will be listed on biowasm.com (optional)
```

## Clean up repos

Tool sources are disposable - `bin/fetch.sh` re-clones whatever the next build needs:

```bash
rm -rf tools/*/src tools/hmmer/easel tools/hmmer/simde
```

## Clean up GitHub Actions

```bash
$ gh run list --limit 200 --repo biowasm/biowasm --json databaseId  -q '.[].databaseId' |
  xargs -IID gh api \
    "repos/$(gh repo view --json nameWithOwner -q .nameWithOwner)/actions/runs/ID" \
    -X DELETE
```
