## wgsim.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/aioli/latest/aioli.js"></script>
<script>
let wgsim = new Aioli("wgsim/2011.10.17");

document.write("Loading...");
wgsim
    // Initialize wgsim
    .init()
    // Create mock reference fasta file
    .then(() => wgsim.fs("writeFile", "/ref.fa", ">chr1\nACGTACGACTAGCAG\n>chr2\nACGATCATACCAGCA"))
    // Run "wgsim" command (overridding default options since reference is very small)
    .then(() => wgsim.exec("-N 5 -d 2 -s 1 -1 5 -2 5 /ref.fa /read1.fastq /read2.fastq"))
    // Get the file contents of one of the two FASTQ files generated
    .then(() => wgsim.cat("/read2.fastq"))
    // Output contents to screen
    .then(d => document.write(`<pre>${d}</pre>`));
</script>
```

### Patch
- Removes a few lines of code that generate unused stdout, which introduces very expensive JS <--> Wasm calls
