![biowasm logo](https://avatars.githubusercontent.com/u/62475458?s=200&v=4)

# biowasm

![Tests](https://github.com/biowasm/biowasm/workflows/Tests/badge.svg)
![Deploy](https://github.com/biowasm/biowasm/workflows/Deploy%20biowasm%20v3/badge.svg)

A repository of genomics tools, compiled from C, C++ and Rust to WebAssembly so they can run in a web browser.

## Getting started

* [Documentation](https://biowasm.com/documentation)
* [List of supported packages](https://biowasm.com/cdn)

## Who uses biowasm?    

| Tool | Why biowasm? |
|-|-|
| [sandbox.bio](https://sandbox.bio) | Runs command-line tools in the browser to power interactive tutorials |
| [42basepairs](https://42basepairs.com) | Runs `samtools`, `bedtools`, `bcftools` and other tools to preview genomic files |
| [SeqtoID](https://seqtoid.org/) ([repo](https://github.com/IT-Academic-Research-Services/seqtoid-web)) | Runs `htsfile` and `seqtk` to identify data issues before file upload |
| [Nanopore](https://bonito.epi2me.io) | In-browser basecaller, uses `minimap2` and `samtools` after basecalling |
| [ViralWasm](https://niema-lab.github.io/ViralWasm/) ([repo](https://github.com/Niema-Lab/ViralWasm)) | Runs `minimap2` and `ViralConsensus` for viral molecular epidemiology analysis |
| [CRISPRLungo](https://pinellolab.github.io/CRISPRLungo/) ([repo](https://github.com/pinellolab/CRISPRLungo)) | Runs `minimap2` as part of a long-read CRISPR analysis pipeline |
| [Datagrok](https://datagrok.ai) ([repo](https://github.com/datagrok-ai/public)) | Runs `kalign` in the browser for multiple-sequence alignment analysis |
| [bedqc](https://quinlan-lab.github.io/bedqc) ([repo](https://github.com/quinlan-lab/bedqc)) | Runs `bedtools` in the browser to validate BED files |
| [Ribbon](https://genomeribbon.com) ([repo](https://github.com/MariaNattestad/Ribbon)) | Runs `samtools` in the browser to parse, estimate coverage and subsample BAM files |
| [fastq.bio](https://www.fastq.bio) ([repo](https://github.com/robertaboukhalil/fastq.bio)) | Runs `fastp` in the browser to evaluate sequencing data quality |
| [GBIF](https://www.gbif.org/occurrence/search?view=dashboard&nucleotideSequence.sequence=TTATATTTTATATTCGGAATTTGGGCAGGAATAATTGGGACATCATTAAGATTATTAATTCGAGCTGAATTAGGAAACCCAGGATCTCTAATTGGGGATGATCAAATTTATAACTCAATTGTAACAGCTCATGCTTTTATTATAATTTTCTTTATAGTAATACCTATTATAATTGGAG&nucleotideSequence.similarity=%5B%2298%22%2C%2299%22%5D) ([repo](https://github.com/gbif/gbif-web)) | Runs `kalign` to align sequences for its sequence phylogeny dashboard charts |
| [Datamonkey](https://v3.datamonkey.org) ([repo](https://github.com/veg/datamonkey3)) | Runs `hyphy` in the browser for molecular evolution analysis |
| [json4u](https://json4u.com) ([repo](https://github.com/loggerhead/json4u)) | Runs `jq` in the browser to power its JSON query feature |
| [Sequence Manipulation Suite](https://paulstothard.github.io/sequence-manipulation-suite/) ([repo](https://github.com/paulstothard/sequence-manipulation-suite)) | Runs `minimap2`, `muscle`, `samtools` and other tools in the browser alongside its classic sequence analysis utilities |
| [Wasm4Learn](https://ifb-elixirfr.github.io/Wasm4Learn/) ([repo](https://github.com/IFB-ElixirFr/Wasm4Learn-content)) | Runs `gawk`, `grep` and `sed` in the browser for bioinformatics teaching |
| [Diana](https://diana.cpg.org.au) ([repo](https://github.com/Wytamma/diana)) | Runs `samtools`, `bgzip` and `tabix` in the browser to index genomes and annotations for TraDIS analysis |

## How it works

| Tool | Description | Link |
|-|-|-|
| biowasm | Recipes for compiling C, C++ and Rust genomics tools to WebAssembly | This repo |
| biowasm CDN | Free server hosting pre-compiled tools for use in your apps | [biowasm.com/cdn](https://biowasm.com/cdn) |
| Aioli | Tool for running these modules in a browser, inside WebWorkers | [biowasm/aioli](https://github.com/biowasm/aioli) |


## Logo

* Logo by [tinygraphs](https://www.tinygraphs.com/#?name=biowasm&shape=labs%2Fisogrids%2Fhexa&theme=seascape&numcolors=4#tryitout)

## Contributing

See [CONTRIBUTING.md](https://github.com/biowasm/biowasm/blob/main/CONTRIBUTING.md).
