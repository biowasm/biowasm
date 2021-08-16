## bhtsne.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
let tsne = new Aioli("bhtsne/2016.08.22");

let output = await CLI.exec(`bhtsne /bhtsne/brain8.snd`)
document.write(`<pre>${output}</pre>`);
</script>
```
