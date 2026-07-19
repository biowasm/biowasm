mkdir src/example
cp test_data/pseudocat.fa test_data/pseudopig2.fa src/example/

# clang 23 (Emscripten 6) adds new warnings (e.g. -Wunused-but-set-global on lastz's
# showProgress/dbgInhibitSegmentReduction) that lastz's -Werror promotes to hard errors.
# Drop -Werror so upstream's harmless warnings don't fail the build.
sed -i 's/-Werror//g' src/Makefile

emmake make CC=emcc LDFLAGS="$EM_FLAGS -O3 --preload-file example@/lastz/" build_lastz
mv src/{lastz,lastz_D}.{mjs,wasm,data} ../build/
