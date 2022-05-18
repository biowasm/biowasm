## htslib.wasm

Used by samtools.

### Usage

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
// There are 3 utils included in htslib:
const CLI = await new Aioli([
	"htslib/tabix/1.10",
	"htslib/htsfile/1.10",
	"htslib/bgzip/1.10"
]);

let output = "";
output += await CLI.exec(`tabix --help`);
output += await CLI.exec(`htsfile --help`);
output += await CLI.exec(`bgzip --help`);

document.write(`<pre>${output}</pre>`);
</script>
```

### Patch
- Making sure the output is a `.js` in the `build` folder
