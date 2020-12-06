## bowtie2.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/aioli/latest/aioli.js"></script>
<script>
let bowtie2 = new Aioli("bowtie2/2.4.2");

document.write("Loading...");
bowtie2
    // Initialize bowtie2
    .init()
    // Run "bowtie" command to map reads to the lambda genome
    .then(() => bowtie2.exec("-x /bowtie2/example/index/lambda_virus -U /bowtie2/example/reads/reads_1.fq"))
    // Output result
    .then(d => document.write(`<pre>${d.stdout}\n${d.stderr}</pre>`));
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
