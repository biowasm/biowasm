## bowtie2.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
let CLI = await new Aioli("bowtie2/bowtie2-align-s/2.4.2", {
  printInterleaved: false
});

// Map reads to the lambda genome
let ref = "/bowtie2/example/index/lambda_virus";
let reads = "/bowtie2/example/reads/reads_1.fq";
let output = await CLI.exec(`bowtie2-align-s -x ${ref} -U ${reads}`)

document.write(`<pre>${output.stdout}</pre>`);
</script>
```


### Patch
- Disable SIMD (simulates it with SIMDe)
- Disable memory-mapping for now
- Disable threads + call `multiseedSearchWorker( (void*)&tids.back() );` manually to replace `threads.push_back()`

### To do
- Add thread support
- Test effect of memory-mapped files
- Write `sse_wrap.h` equivalent for WebAssembly's SIMD instructions?
