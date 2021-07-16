## ssw.wasm

### Warning

This implementation of Smith-Waterman relies on SIMD instructions to run in a performant way. If the end user does not have SIMD enabled, Aioli will load the `ssw-nosimd.wasm` version, though it is considerably slower.

To enable SIMD in your browser:

* **Firefox**: Go to [about:config](about:config), search for `javascript.options.wasm_simd` and toggle to `true` (enabled by default as of Firefox 89)
* **Chrome**: Go to [chrome://flags](chrome://flags), search for `WebAssembly SIMD support` and toggle to `Enabled` (enabled by default as of Chrome 91)

### Usage

```html
<script src="https://cdn.biowasm.com/aioli/latest/aioli.js"></script>
<script>
let ssw = new Aioli("ssw/1.2.4");

document.write("Loading...");
ssw
    // Initialize ssw
    .init()
    // Create mock reference & query fasta files
    .then(() => ssw.fs("writeFile", "/ref.fa", ">chr1\nACTACGACTACGACTACGACGCGCGATTCGCGCGCCGATATACGACTACGACTA\n"))
    .then(() => ssw.fs("writeFile", "/query.fa", ">read1\nTACGACTACG\n"))
    // Run "ssw" command to align query to reference
    .then(() => ssw.exec("-h -s -c ref.fa query.fa"))
    // Output resulting SAM file
    .then(d => document.write(`<pre>${d.stdout}</pre>`));
</script>
```
