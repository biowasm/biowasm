![biowasm logo](https://avatars.githubusercontent.com/u/62475458?s=200&v=4)

# biowasm

![cdn-stg.biowasm.com](https://github.com/biowasm/biowasm/workflows/Deploy%20biowasm%20v2/badge.svg)

A repository of genomics tools, compiled from C/C++ to WebAssembly so they can run in a web browser.

## Getting started

1. Visit [biowasm.com](https://biowasm.com/) to see how you can start using biowasm with just a few lines of code.
2. Check out our [Getting Started](https://github.com/biowasm/aioli#a-simple-example) guide.

## Who uses biowasm?    

| Tool | Why biowasm? |
|-|-|
| [sandbox.bio](https://sandbox.bio) | Run command-line tools in the browser to power interactive tutorials |
| [bedqc](https://quinlan-lab.github.io/bedqc) ([repo](https://github.com/quinlan-lab/bedqc)) | Run `bedtools` in the browser to validate BED files |
| [Ribbon](https://genomeribbon.com) ([repo](https://github.com/MariaNattestad/Ribbon)) | Run `samtools` in the browser to parse, estimate coverage and subsample BAM files |
| [fastq.bio](https://www.fastq.bio) ([repo](https://github.com/robertaboukhalil/fastq.bio)) | Run `fastp` in the browser to evaluate sequencing data quality |

## Supported tools

C/C++ tools that have been compiled to WebAssembly:

#### Bioinformatics tools

| Tool | Version | Description |
|-|-|-|
| [bedtools](tools/bedtools2) | 2.29 | Parse `.bed` files and perform complex "genome arithmetic" |
| [samtools](tools/samtools) | 1.10 | Parse `.sam` / `.bam` read alignment files |
| [bcftools](tools/bcftools) | 1.10 | Parse `.vcf` / `.bcf` variant calling files |
| [htslib](tools/htslib) | 1.10 | Bioinformatics file format utilities: `tabix`, `htsfile`, and `bgzip` |
| [seqtk](tools/seqtk) | 1.3 | Parse and wrangle `.fasta` / `.fastq` files |
| [fastp](tools/fastp) | 0.20.1 | Evaluate data quality of `.fastq` files |
| [bowtie2](tools/bowtie2) | 2.4.2 | Align sequencing reads (`.fastq`) files to a reference genome |
| [minimap2](tools/minimap2) | 2.22 | Pairwise sequence alignment |
| [kalign](tools/kalign) | 3.3.1 | Multiple sequence alignment |
| [fasttree](tools/fasttree) | 2.1.11 | Build phylogenetic trees from multiple sequence alignments |
| [seq-align](tools/seq-align) | 2017.10.18 | Align sequences using Smith-Waterman/Needleman-Wunsch |
| [ssw](tools/ssw) | 1.2.4 | A SIMD implementation of the Smith-Waterman algorithm |
| [ivar](tools/ivar) | 1.3.1 | Tools for viral amplicon-based sequencing |
| [modbam2bed](tools/modbam2bed) | 0.3.1 | Summarize `.bam` modified bases as `.bed` files with counts (thanks [@cjw85](https://github.com/cjw85)) |
| [wgsim](tools/wgsim) | 2011.10.17 | Simulate short reads from a reference genome |
| [bhtsne](tools/bhtsne) | 2016.08.22 | Run the t-SNE dimensionality-reduction algorithm |

#### General utilities

| Tool | Version | Description |
|-|-|-|
| [jq](tools/jq) | 1.6 | Filter and wrangle `JSON` strings |
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

## Contributing

See [CONTRIBUTING.md](https://github.com/biowasm/biowasm/blob/main/CONTRIBUTING.md).
