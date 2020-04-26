# biowasm
WebAssembly modules for common genomics utilities, including:

* samtools v1.10 (and htslib)
* bedtools v2.29
* fastp v0.20.1
* seqtk v1.3
* wgsim
* bhtsne
* seq-align


## Setup

```bash
# Emscripten version to use (most tools were tested with 1.39.1)
TAG=1.39.1

# Fetch Emscripten docker image
docker pull robertaboukhalil/emsdk:$TAG

# Create the container and mount ~/wasm to /src in the container
docker run \
    -dt \
    -p 12345:80 \
    --name wasm \
    --volume ~/wasm:/src \
    robertaboukhalil/emsdk:$TAG
```


## Compile a tool

```bash
# Go into your container
docker exec -it wasm bash

# Compile seqtk
cd biowasm/
make seqtk

# This will create tools/<tool name>/build with .js/.wasm files
ls tools/seqtk/build
```


## Contribute new tool

First, add the tool as a git module:

```bash
# Fetch codebase
mkdir -p tools/seqtk
git submodule add https://github.com/lh3/seqtk.git tools/seqtk/src

# Get specific version of the tool
cd tools/seqtk/src
git checkout v1.3
cd -

# Stage changes for git
git add tools/seqtk/src .gitmodules
```

You should also create the following files:

```bash
tools/<tool>/
    README.md   Details about the tool and dependencies
    compile.sh  Script that will run to compile the tool to WebAssembly (can use `$EM_FLAGS` for common flags)
    patch       Patches that need to be applied to the code to compile it to WebAssembly (optional)
```

## Todo

- Add tests
- Add support for compiling bioinformatics tools written in Rust such as [sourmash](https://github.com/dib-lab/sourmash/tree/v3.2.2/src/core) and [rust-bio](https://github.com/rust-bio/rust-bio)

