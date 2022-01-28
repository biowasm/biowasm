## grep.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
// Note that we use `script type="module"` so we can use top-level await statements
const CLI = await new Aioli("grep/3.7");

// Create sample file
await CLI.fs.writeFile("data.txt", "hello\ngood\nmorning\n");

// Retrieve lines that start with g
let output = await CLI.exec("grep ^g data.txt");
document.write(`<pre>${output}</pre>`);
</script>
```
