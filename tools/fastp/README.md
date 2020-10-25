## fastp.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/aioli/latest/aioli.js"></script>
<script>
let fastp = new Aioli("fastp/0.20.1");

document.write("Loading...");
fastp
    // Initialize fastp
    .init()
    // Run "fastp" command on a sample FASTQ file
    .then(() => fastp.exec("-i /fastp/testdata/R1.fq"))
    // Output resulting JSON file that contains metrics
    .then(() => fastp.cat("/fastp.json"))
    .then(d => document.write(`<pre>${d}</pre>`));
</script>
```

### Patch
- Disables threads to make it easier to compile to WebAssembly (need to update Aioli to support Wasm modules with threads)

### TODO
- Enable threads
