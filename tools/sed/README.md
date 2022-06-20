## sed.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
const CLI = await new Aioli("sed/4.8");

// Create mock data
await CLI.fs.writeFile("test.fastq", "@read1\nACGTACGACTAGCAG\n+\nJJJJJJJJJJJJJJJ\n@read2\nACGATCATACCAGCA\n+\nJJJJJJJJJJJJJJJ\n");

// Simple find/replace
const output = await CLI.exec("sed s/GACT/----/ test.fastq");
document.write(`<h3>Find/replace</h3><pre>${output}</pre>`);

// Convert FASTQ to FASTA (from https://github.com/stephenturner/oneliners#awk--sed-for-bioinformatics)
const fasta = await CLI.exec("sed -n 1~4s/^@/>/p;2~4p test.fastq");
document.write(`<h3>FASTQ to FASTA</h3><pre>${fasta}</pre>`);
</script>
```
