### Emscripten
- `docker pull robertaboukhalil/emsdk:1.39.1`

### Dependencies
- Run `make htslib` before running `make samtools`

### Patch
- Need to reset `opt` variables so that it works properly when call `main()` multiple times
