## wgsim.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
let CLI = await new Aioli("wgsim/2011.10.17");

// Create mock reference fasta file
await CLI.fs.writeFile("ref.fa", ">chr1\nACGTACGACTAGCAG\n>chr2\nACGATCATACCAGCA");

// Run "wgsim" command (overridding default options since reference is very small)
// and output result to read1.fastq and read2.fastq
await CLI.exec("wgsim -N 5 -d 2 -s 1 -1 5 -2 5 ref.fa read1.fastq read2.fastq");

// Get the file contents of one of the two FASTQ files generated
const output = await CLI.cat("read2.fastq");
document.write(`<pre>${output}</pre>`);
</script>
```

### Patch
- Removes a few lines of code that generate unused stdout, which introduces very expensive JS <--> Wasm calls
