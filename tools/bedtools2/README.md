## bedtools2.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/aioli/latest/aioli.js"></script>
<script>
let bedtools2 = new Aioli("bedtools2/2.29.2");

document.write("Loading...");
bedtools2
    // Initialize bedtools2
    .init()
    // Run "bedtools" command to intersect 2 BED files
    .then(() => bedtools2.exec("intersect -a /bedtools2/test/intersect/a.bed -b /bedtools2/test/intersect/b.bed"))
    // Output result
    .then(d => document.write(`<pre>${d.stdout}\n${d.stderr}</pre>`));
</script>
```


### Patch
- Use em++ instead of g++
