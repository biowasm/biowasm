![biowasm logo](https://avatars.githubusercontent.com/u/62475458?s=200&v=4)

# biowasm

![Tests](https://github.com/biowasm/biowasm/workflows/Tests/badge.svg)
![Deploy](https://github.com/biowasm/biowasm/workflows/Deploy%20biowasm%20v3/badge.svg)

A repository of genomics tools, compiled from C/C++ to WebAssembly so they can run in a web browser.

## Getting started

* [Documentation](https://biowasm.com/documentation)
* [List of supported packages](https://biowasm.com/cdn)

## Who uses biowasm?    

| Tool | Why biowasm? |
|-|-|
| [sandbox.bio](https://sandbox.bio) | Run command-line tools in the browser to power interactive tutorials |
| [bedqc](https://quinlan-lab.github.io/bedqc) ([repo](https://github.com/quinlan-lab/bedqc)) | Run `bedtools` in the browser to validate BED files |
| [Ribbon](https://genomeribbon.com) ([repo](https://github.com/MariaNattestad/Ribbon)) | Run `samtools` in the browser to parse, estimate coverage and subsample BAM files |
| [fastq.bio](https://www.fastq.bio) ([repo](https://github.com/robertaboukhalil/fastq.bio)) | Run `fastp` in the browser to evaluate sequencing data quality |

## How it works

| Tool | Description | Link |
|-|-|-|
| biowasm | Recipes for compiling C/C++ genomics tools to WebAssembly | This repo |
| biowasm CDN | Free server hosting pre-compiled tools for use in your apps | [biowasm.com/cdn](https://biowasm.com/cdn) |
| Aioli | Tool for running these modules in a browser, inside WebWorkers | [biowasm/aioli](https://github.com/biowasm/aioli) |


## Logo

* Logo by [tinygraphs](https://www.tinygraphs.com/#?name=biowasm&shape=labs%2Fisogrids%2Fhexa&theme=seascape&numcolors=4#tryitout)

## Contributing

See [CONTRIBUTING.md](https://github.com/biowasm/biowasm/blob/main/CONTRIBUTING.md).
