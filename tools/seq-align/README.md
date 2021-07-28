## seq-align.wasm

### Usage

Note that `seq-align` has multiple sub-tools so you need to specify which one you want Aioli to load in the constructor:

#### Smith-Waterman

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
let CLI = await new Aioli("seq-align/smith_waterman/2017.10.18");

// Align two sequences
let output = await CLI.exec("smith_waterman ACGT ACCCCGT");
document.write(`<pre>${output}</pre>`);
</script>
```

#### Needleman-Wunsch

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
let CLI = await new Aioli("seq-align/needleman_wunsch/2017.10.18");

// Align two sequences
let output = await CLI.exec("needleman_wunsch ACGT ACCCCGT");
document.write(`<pre>${output}</pre>`);
</script>
```

### Patch
- Add `-s USE_ZLIB=1` for gzipped file support
- Replace `ar` with `$(AR)` to make sure we pick up the Emscripten ar
