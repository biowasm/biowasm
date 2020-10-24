## seqtk.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/aioli/latest/aioli.js"></script>
<script>
let seqtk = new Aioli("seqtk/1.3");

document.write("Loading...");
seqtk
    // Initialize seqtk
    .init()
    // Create mock fasta file
    .then(() => seqtk.fs("writeFile", "/test.fa", ">chr1\nACGTACGACTAGCAG\n>chr2\nACGATCATACCAGCA"))
    // Run "seqtk comp" command to calculate nucleotide composition
    .then(() => seqtk.exec("comp /test.fa"))
    // Output result
    .then(d => document.write(`<pre>${d.stdout}\n${d.stderr}</pre>`));
</script>
