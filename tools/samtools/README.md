## samtools.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/aioli/latest/aioli.js"></script>
<script>
let samtools = new Aioli("samtools/1.10");

document.write("Loading...");
samtools
    // Initialize samtools
    .init()
    // Run "samtools view" command with "-q 20" filter
    .then(() => samtools.exec("view -q 20 /samtools/examples/toy.sam"))
    // Output result
    .then(d => document.write(`<pre>${d.stdout}\n${d.stderr}</pre>`));
</script>
```

### Dependencies
- Run `make htslib` before running `make samtools`

### Patch
- Need to reset `opt` variables so that it works properly when call `main()` multiple times
