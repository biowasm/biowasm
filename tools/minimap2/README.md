## minimap2.wasm

### Warning

Minimap2 relies on SIMD instructions. If the end user does not have SIMD enabled, Aioli will load the version without SIMD, though it is considerably slower.

To enable SIMD in your browser:

* **Firefox**: Go to [about:config](about:config), search for `javascript.options.wasm_simd` and toggle to `true` (enabled by default as of Firefox 89)
* **Chrome**: Go to [chrome://flags](chrome://flags), search for `WebAssembly SIMD support` and toggle to `Enabled` (enabled by default as of Chrome 91)

### Usage

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
let CLI = await new Aioli("minimap2/2.22");

const output = await CLI.exec("minimap2 -a /minimap2/MT-human.fa /minimap2/MT-orang.fa");
document.write(`<pre>${output}</pre>`);
</script>
```
