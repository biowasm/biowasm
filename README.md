# biowasm

![cdn-stg.biowasm.com](https://github.com/biowasm/biowasm/workflows/Deploy%20biowasm-stg/badge.svg) ![cdn.biowasm.com](https://github.com/biowasm/biowasm/workflows/Deploy%20biowasm-prd/badge.svg)

A repository of genomics tools, compiled from C/C++ to WebAssembly so they can run in a web browser.

## Getting started

Check out our [Getting Started](https://github.com/biowasm/aioli#a-simple-example) guide.

## Supported tools

C/C++ tools that have been compiled to WebAssembly:

| Tool | Version | Description |
|-|-|-|
| [samtools/htslib](tools/samtools) | 1.10 | Parse and manipulate <code>.sam</code> / <code>.bam</code> read alignment files |
| [bedtools2](tools/bedtools2) | 2.29 | Parse <code>.bed</code> files and perform complex "genome arithmetic" |
| [bowtie2](tools/bowtie2) | 2.4.2 | Align sequencing reads (<code>.fastq</code>) files to a reference genome |
| [fastp](tools/fastp) | 0.20.1 | Manipulate and evaluate QC of <code>.fastq</code> files |
| [seqtk](tools/seqtk) | 1.3 | Manipulate and evaluate QC of <code>.fasta</code> / <code>.fastq</code> files |
| [ssw](tools/ssw) | 1.2.4 | A SIMD implementation of the Smith-Waterman algorithm |
| [wgsim](tools/wgsim) | 2011.10.17 | Simulate short reads from a reference genome |
| [seq-align](tools/seq-align) | 2017.10.18 | Align sequences using Smith-Waterman/Needleman-Wunsch algorithms |
| [bhtsne](tools/bhtsne) | 2016.08.22 | Run the t-SNE dimensionality-reduction algorithm |

## How it works

| Tool | Description | Link |
|-|-|-|
| biowasm | Recipes for compiling C/C++ genomics tools to WebAssembly | This repo |
| biowasm CDN | Free CDN server hosting pre-compiled tools for use in your apps | See [cdn.biowasm.com](https://cdn.biowasm.com) |
| Aioli | Tool for running these modules in a browser, inside WebWorkers | See [biowasm/aioli](https://github.com/biowasm/aioli) |


## Tools using biowasm

| Tool | URL | Repo |
|-|-|-|
| Ribbon | [genomeribbon.com](https://genomeribbon.com) | [MariaNattestad/Ribbon](https://github.com/MariaNattestad/Ribbon) |
| Alignment Sandbox | [alignment.sandbox.bio](https://alignment.sandbox.bio/) | [RobertAboukhalil/alignment-sandbox](https://github.com/robertaboukhalil/alignment-sandbox) |
| tSNE Sandbox | [tsne.sandbox.bio](https://tsne.sandbox.bio/) | [RobertAboukhalil/tsne-sandbox](https://github.com/robertaboukhalil/tsne-sandbox) |
| fastq.bio | [fastq.bio](http://www.fastq.bio/) | [RobertAboukhalil/fastq.bio](https://github.com/robertaboukhalil/fastq.bio) |
| bam.bio | [bam.bio](http://www.bam.bio/) | [RobertAboukhalil/bam.bio](https://github.com/robertaboukhalil/bam.bio) |

---

## Contributing

Ignore the rest of this README if you are not contributing changes to the biowasm repo.

### Setup

Tools listed in biowasm were compiled to WebAssembly using `Emscripten 2.0.25`.

```bash
# Fetch Emscripten docker image
docker pull emscripten/emsdk:2.0.25

# Create the container and mount ~/wasm to /src in the container
docker run \
    -it -d \
    -p 80:80 \
    --name wasm \
    --volume ~/wasm:/src \
    emscripten/emsdk:2.0.25

# Go into the container
docker exec -u root -it wasm bash
# While inside the container, install dependencies
apt-get update
apt-get install -y autoconf liblzma-dev less vim
# Create small web server for testing
cat << EOF > server.py
import http.server
import socketserver

handler = http.server.SimpleHTTPRequestHandler
handler.extensions_map['.wasm'] = 'application/wasm'
httpd = socketserver.TCPServer(('', 80), handler)
httpd.serve_forever()
EOF
chmod +x server.py
# Launch the web server
python3.7 /src/server.py &
```


### Compile a tool

```bash
# Go into your container
docker exec -it wasm bash

# Set up biowasm (only need to do this once)
cd biowasm/
make init

# Compile seqtk
VERSION=1.2 BRANCH=v1.2 make seqtk

# This will create tools/<tool name>/build with .js/.wasm files
ls tools/seqtk/build
```


### Add a new tool

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
    README.md        Details about the tool and dependencies
    compile.sh       Script that will run to compile the tool to WebAssembly (can use `$EM_FLAGS` for common flags)
    patches/    
        <tag>        Patch applied to the code to compile it to WebAssembly; branch- or tag-specific (optional)
    configs/
        <tag>.json   Configuration file with info about which WebAssembly features are needed (see ssw for an example); branch- or tag-specific (optional)
```

## Deploy changes

* Changes merged are auto-deployed via GitHub Actions to `cdn-stg.biowasm.com/v2`.


## To do

- Deploy one tool without re-compiling all others: download data from the CDN onto the GitHub Actions VM first?
- Run each tool's tests: use Selenium? Can't use node.js when have `.data` files
- Generate HTML file for each tool: CLI for testing, predefined queries, etc
- Support for Rust bioinformatics tools such as [sourmash](https://github.com/dib-lab/sourmash/tree/v3.2.2/src/core) and [rust-bio](https://github.com/rust-bio/rust-bio)
