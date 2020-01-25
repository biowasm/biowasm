### Emscripten
- `docker pull robertaboukhalil/emsdk:1.39.1`

### Patch
- Add `-s USE_ZLIB=1` for gzipped file support
- Replace `ar` with `$(AR)` to make sure we pick up the Emscripten ar
