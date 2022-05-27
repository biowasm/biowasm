## ivar.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
let CLI = await new Aioli(["samtools/1.10", "ivar/1.3.1"]);

// Create a mock BAM file
await CLI.exec(`samtools view -bS -o test.bam /samtools/examples/toy.sam`);

// Use ivar to trim primers
await CLI.exec(`ivar trim -i test.bam -p test.trimmed -m 5`);

// Output trimmed BAM file
let output = await CLI.exec(`samtools view test.trimmed.bam`);
document.write(`<pre>${output}</pre>`);
</script>
```

### Patch
- Reset `optind` and CLI params so that we can call `main()` multiple times in a row
