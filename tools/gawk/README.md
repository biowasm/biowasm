## gawk.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
// Note that we use `script type="module"` so we can use top-level await statements
const CLI = await new Aioli("gawk/5.1.0");

// Create mock tab separated file
await CLI.fs.writeFile("data.tsv", `column1\tcolumn2\tcolumn3\n1\t2\t3\n4\t5\t6\n7\t8\t9\n`);

// Retrieve 2nd column
let output = await CLI.exec("gawk", [ 'BEGIN{ sum = 0; }{ print $2; sum += $2; } END { print("Total of column 2 = " sum) }', "data.tsv" ]);
document.write(`<pre>${output}</pre>`);
</script>
```

### Patch

- Avoid printing signal warnings to stderr
