## bhtsne.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/aioli/latest/aioli.js"></script>
<script>
let tsne = new Aioli("bhtsne/t-sne/2016.08.22");

document.write("Loading (may take a few seconds)...");
tsne
    // Initialize tsne
    .init()
    // Run dimensionality reduction algorithm on example matrix 
    .then(() => tsne.exec("/bhtsne/brain8.snd"))
    // Output result
    .then(d => document.write(`<pre>${d.stdout}\n${d.stderr}</pre>`));
</script>
```
