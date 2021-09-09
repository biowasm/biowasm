## jq.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
let CLI = await new Aioli("jq/1.6");

// Create mock JSON file
await CLI.fs.writeFile("test.json", `[{"some":{"data": "here"}},{"some": {"data": "there"}}]`);

// Retrieve JSON subset
let output = await CLI.exec("jq", [
	"-r", ".[1].some.data", "test.json"
]);
document.write(`<pre>${output}</pre>`);
</script>
```

### Patch
- Reset `options` variable so that we can call `main()` repeatedly
