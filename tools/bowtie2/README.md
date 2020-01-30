### Emscripten
- `docker pull robertaboukhalil/emsdk:1.39.1`

### Patch
- Disable SIMD (simulates it with SIMDe)
- Disable memory-mapping for now
- Disable threads + call `multiseedSearchWorker( (void*)&tids.back() );` manually to replace `threads.push_back()`

### TODO
- Add thread support
- Test effect of memory-mapped files
- Write `sse_wrap.h` equivalent for WebAssembly's SIMD instructions?
