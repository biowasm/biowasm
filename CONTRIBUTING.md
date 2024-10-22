# Contributing

## Development Setup

Tools listed in biowasm were compiled to WebAssembly using `Emscripten 2.0.25`. Here is how to set up your dev environment:

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

First, add the tool as a git module:

```bash
TOOL=minimap2
BRANCH=v2.22
REPO=https://github.com/lh3/minimap2.git

# Fetch codebase
mkdir -p tools/$TOOL
git submodule add $REPO tools/$TOOL/src

# Get specific version of the tool
cd tools/$TOOL/src
git checkout $BRANCH
git submodule update --init --recursive
cd -

# Stage changes for git
git add tools/$TOOL/src .gitmodules

# Then create a compile.sh script that compiles the tool to WebAssembly
echo "# TODO" > tools/$TOOL/compile.sh
chmod +x tools/$TOOL/compile.sh

mkdir -p tools/$TOOL/examples
echo "#TODO" > tools/$TOOL/examples/$BRANCH.html
```

You should also modify/create the following files:

```bash
biowasm.json         Add the new tool here so it gets deployed to the CDN

tools/<tool>/
    compile.sh       Script that will run to compile the tool to WebAssembly (use `$EM_FLAGS` for common flags)
    patches/    
        <tag>.patch  Patch applied to the code to compile it to WebAssembly; branch- or tag-specific (optional)
    examples/
        <tag>.html   Example HTML code for running the tool; will be listed on biowasm.com (optional)
```

## Clean up GitHub Actions

```bash
$ gh run list --limit 200 --repo biowasm/biowasm --json databaseId  -q '.[].databaseId' |
  xargs -IID gh api \
    "repos/$(gh repo view --json nameWithOwner -q .nameWithOwner)/actions/runs/ID" \
    -X DELETE
```
