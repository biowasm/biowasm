# biowasm

![cdn-stg.biowasm.com](https://github.com/biowasm/biowasm/workflows/Deploy%20biowasm%20v2/badge.svg)

A repository of genomics tools, compiled from C/C++ to WebAssembly so they can run in a web browser.

## Getting started

1. Visit [biowasm.com](https://biowasm.com/) to see how you can start using biowasm with just a few lines of code.
2. Check out our [Getting Started](https://github.com/biowasm/aioli#a-simple-example) guide.

## Who uses biowasm?    

| Tool | URL | Repo |
|-|-|-|
| sandbox.bio | [sandbox.bio](https://sandbox.bio) | - |
| 42basepairs | [42basepairs.com](https://42basepairs.com) | - |
| Ribbon | [genomeribbon.com](https://genomeribbon.com) | [MariaNattestad/Ribbon](https://github.com/MariaNattestad/Ribbon) |
| BEDQC | [quinlan-lab.github.io/bedqc](https://quinlan-lab.github.io/bedqc/) | [Quinlan-Lab/BEDQC](https://github.com/quinlan-lab/bedqc) |
| fastq.bio | [fastq.bio](http://www.fastq.bio) | [RobertAboukhalil/fastq.bio](https://github.com/robertaboukhalil/fastq.bio) |
| tSNE Sandbox | [tsne.sandbox.bio](https://tsne.sandbox.bio) | [RobertAboukhalil/tsne-sandbox](https://github.com/robertaboukhalil/tsne-sandbox) |
| Alignment Sandbox | [alignment.sandbox.bio](https://alignment.sandbox.bio) | [RobertAboukhalil/alignment-sandbox](https://github.com/robertaboukhalil/alignment-sandbox) |

## Supported tools

C/C++ tools that have been compiled to WebAssembly:

#### Bioinformatics tools

| Tool | Version | Description |
|-|-|-|
| [samtools](tools/samtools) | 1.10 | Parse and manipulate <code>.sam</code> / <code>.bam</code> read alignment files |
| [bcftools](tools/bcftools) | 1.10 | Parse and manipulate <code>.vcf</code> / <code>.bcf</code> variant calling files |
| [bedtools](tools/bedtools2) | 2.29 | Parse <code>.bed</code> files and perform complex "genome arithmetic" |
| [bowtie2](tools/bowtie2) | 2.4.2 | Align sequencing reads (<code>.fastq</code>) files to a reference genome |
| [minimap2](tools/minimap2) | 2.22 | Align sequences to each other |
| [kalign](tools/kalign) | 3.3.1 | Multiple sequence alignment |
| [fastp](tools/fastp) | 0.20.1 | Manipulate and evaluate QC of <code>.fastq</code> files |
| [seqtk](tools/seqtk) | 1.3 | Manipulate and evaluate QC of <code>.fasta</code> / <code>.fastq</code> files |
| [ssw](tools/ssw) | 1.2.4 | A SIMD implementation of the Smith-Waterman algorithm |
| [modbam2bed](tools/modbam2bed) | 0.3.1 | Summarize <code>.bam</code> files with modified bases as <code>.bed</code> files with counts |
| [wgsim](tools/wgsim) | 2011.10.17 | Simulate short reads from a reference genome |
| [seq-align](tools/seq-align) | 2017.10.18 | Align sequences using Smith-Waterman/Needleman-Wunsch algorithms |
| [bhtsne](tools/bhtsne) | 2016.08.22 | Run the t-SNE dimensionality-reduction algorithm |

#### General utilities

| Tool | Version | Description |
|-|-|-|
| [jq](tools/jq) | 1.6 | Filter and wrangle <code>JSON</code> strings |
| [gawk](tools/gawk) | 5.1.0 | Manipulate data files with patterns of interest |
| [grep](tools/grep) | 3.7 | Search and filter files |


## How it works

| Tool | Description | Link |
|-|-|-|
| biowasm | Recipes for compiling C/C++ genomics tools to WebAssembly | This repo |
| biowasm CDN | Free server hosting pre-compiled tools for use in your apps | [cdn.biowasm.com](https://cdn.biowasm.com/v2/) |
| Aioli | Tool for running these modules in a browser, inside WebWorkers | [biowasm/aioli](https://github.com/biowasm/aioli) |


## Logo

* Logo by [tinygraphs](https://www.tinygraphs.com/#?name=biowasm&shape=labs%2Fisogrids%2Fhexa&theme=seascape&numcolors=4#tryitout)

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

# Other files needed
echo "TODO" > tools/$TOOL/README.md
echo "# TODO" > tools/$TOOL/compile.sh
chmod +x tools/$TOOL/compile.sh
```

You should also create the following files:

```bash
tools/<tool>/
    README.md        Details about the tool and dependencies
    compile.sh       Script that will run to compile the tool to WebAssembly (can use `$EM_FLAGS` for common flags)
    patches/    
        <tag>.patch  Patch applied to the code to compile it to WebAssembly; branch- or tag-specific (optional)
    configs/
        <tag>.json   Configuration file with info about which WebAssembly features are needed (see ssw for an example); branch- or tag-specific (optional)
```

Finally, you can edit:

* `config/tools.json` to make sure the new tool gets deployed
* `cloudflare/cdn/public/index.html` to list the new tool on the CDN's [packages page](https://cdn.biowasm.com/v2/)


### Future candidates

- [Sniffles](https://github.com/fritzsedlazeck/Sniffles)
- [FreeBayes](https://github.com/freebayes/freebayes)
- [HISAT2](https://github.com/DaehwanKimLab/hisat2)
