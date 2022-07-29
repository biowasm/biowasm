## ssw.wasm

### Warning

This implementation of Smith-Waterman relies on SIMD instructions to run in a performant way. If the end user does not have SIMD enabled, Aioli will load the non-SIMD version, though it is considerably slower.

To enable SIMD in your browser:

* **Firefox**: Go to [about:config](about:config), search for `javascript.options.wasm_simd` and toggle to `true` (enabled by default as of Firefox 89)
* **Chrome**: Go to [chrome://flags](chrome://flags), search for `WebAssembly SIMD support` and toggle to `Enabled` (enabled by default as of Chrome 91)

### Usage

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
let CLI = await new Aioli("ssw/1.2.4", { printInterleaved: false });

// Create mock reference & query fasta files
await CLI.fs.writeFile("ref.fa", ">chr1\nACTACGACTACGACTACGACGCGCGATTCGCGCGCCGATATACGACTACGACTA\n");
await CLI.fs.writeFile("query.fa", ">read1\nTACGACTACG\n");

// Run "ssw" command to align query to reference and output SAM file
const output = await CLI.exec("ssw -h -s -c ref.fa query.fa");
document.write(`
    stdout: <pre>${output.stdout}</pre>
    <br /><br />
    stderr: <pre>${output.stderr}</pre>
`);
</script>
```
