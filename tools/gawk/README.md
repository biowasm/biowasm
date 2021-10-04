## gawk.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
let CLI = await new Aioli("gawk/5.1.0");
</script>
```

### Patch

- Avoid printing signal warnings to stderr
