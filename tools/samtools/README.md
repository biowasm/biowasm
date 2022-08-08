## samtools.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
let CLI = await new Aioli("samtools/1.10");

// Show reads from toy.sam with flag "16".
// Try replacing "-f 16" with "-f 0".
let file = "/samtools/examples/toy.sam";
let output = await CLI.exec(`samtools view -f 16 ${file}`)
document.write(`<pre>${output}</pre>`);
</script>
```

### Patch
- Need to reset `opt` variables so that it works properly when call `main()` multiple times
- Need to remove `pthread_create` in `samtools sort` so it works without thread support
