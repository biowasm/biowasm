## seq-align.wasm

### Usage

Note that `seq-align` has multiple sub-tools so you need to specify which one you want Aioli to load in the constructor:

#### Smith-Waterman

```html
<script src="https://cdn.biowasm.com/aioli/latest/aioli.js"></script>
<script>
let sw = new Aioli("seq-align/smith_waterman/2017.10.18");

document.write("Loading...");
sw
    // Initialize sw
    .init()
    // Align ACGT and ACCCCGT using Smith-Waterman algorithm
    .then(() => sw.exec("ACGT ACCCCGT"))
    // Output result
    .then(d => document.write(`<pre>${d.stdout}\n${d.stderr}</pre>`));
</script>
```

#### Needleman-Wunsch

```html
<script src="https://cdn.biowasm.com/aioli/latest/aioli.js"></script>
<script>
let nw = new Aioli("seq-align/needleman_wunsch/2017.10.18");

document.write("Loading...");
nw
    // Initialize nw
    .init()
    // Align ACGT and ACCCCGT using Needleman-Wunsch algorithm
    .then(() => nw.exec("ACGT ACCCCGT"))
    // Output result
    .then(d => document.write(`<pre>${d.stdout}\n${d.stderr}</pre>`));
</script>
```

### Patch
- Add `-s USE_ZLIB=1` for gzipped file support
- Replace `ar` with `$(AR)` to make sure we pick up the Emscripten ar
