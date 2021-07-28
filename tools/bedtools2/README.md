## bedtools.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
let CLI = await new Aioli("bedtools/2.29.2");

// Intersect two .bed files
let bed1 = "/bedtools/test/intersect/a.bed";
let bed2 = "/bedtools/test/intersect/b.bed";

let output = await CLI.exec(`bedtools intersect -a ${bed1} -b ${bed2}`);
document.write(`<pre>${output}</pre>`);
</script>
```


### Patch
- Use em++ instead of g++
