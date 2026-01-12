mkdir src/example
cp test_data/pseudocat.fa test_data/pseudopig2.fa src/example/
emmake make CC=emcc LDFLAGS="$EM_FLAGS -O3 --preload-file example@/lastz/" build_lastz
mv src/{lastz,lastz_D}.{js,wasm,data} ../build/
