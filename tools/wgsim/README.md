### Emscripten
- `docker pull robertaboukhalil/emsdk:1.39.1`

### Patch
- Removes a few lines of code that generate unused stdout, which introduces very expensive JS <--> Wasm calls
