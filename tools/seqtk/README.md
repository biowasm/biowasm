## seqtk.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
let CLI = await new Aioli("seqtk/1.3");

// Create mock fasta file
await CLI.fs.writeFile("test.fa", ">chr1\nACGTACGACTAGCAG\n>chr2\nACGATCATACCAGCA");

// Run "seqtk comp" command to calculate nucleotide composition
const output = await CLI.exec("seqtk comp test.fa");
document.write(`<pre>${output}</pre>`);
</script>
```
