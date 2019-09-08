# biowasm
WebAssembly modules for common genomics utilities


```
REPOS=("https://github.com/lh3/seqtk.git" "https://github.com/samtools/samtools.git" "https://github.com/samtools/htslib.git" "https://github.com/arq5x/bedtools2.git")

cd tools/
for REPO in ${REPOS[@]}; do
  git submodule add $REPO
done
```
