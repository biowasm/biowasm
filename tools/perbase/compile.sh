#!/bin/bash

set -e

# rust-htslib runs bindgen over htslib's headers, which needs libclang
command -v clang-20 > /dev/null || {
	sudo apt-get update
	sudo apt-get install -y libclang-20-dev
}
EM_RUST
embuilder build zlib

SYSROOT=$(em-config CACHE)/sysroot
export LIBCLANG_PATH=${LIBCLANG_PATH:-/usr/lib/llvm-20/lib}
export BINDGEN_EXTRA_CLANG_ARGS_wasm32_unknown_emscripten="--sysroot=$SYSROOT -fvisibility=default"
# Rust panics compile to Emscripten exceptions, so the catching they need cannot be disabled
export RUSTFLAGS=$(EM_RUST_FLAGS -O2 -s DISABLE_EXCEPTION_CATCHING=0)

cargo build --release --locked --target wasm32-unknown-emscripten -j4 \
	--config "target.wasm32-unknown-emscripten.z.rustc-link-lib=['z']" \
	--config "target.wasm32-unknown-emscripten.z.include='$SYSROOT/include'" \
	--config "target.wasm32-unknown-emscripten.z.root='$SYSROOT'"
cp target/wasm32-unknown-emscripten/release/perbase.{js,wasm} ../build/
