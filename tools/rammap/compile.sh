#!/bin/bash

set -e

EM_RUST

export RUSTFLAGS="-Ctarget-feature=+simd128 $(EM_RUST_FLAGS -sDEFAULT_TO_CXX)"
cargo build --release --locked --target wasm32-unknown-emscripten -p rammap --bin rammap --no-default-features
cp target/wasm32-unknown-emscripten/release/rammap.{js,wasm} ../build/
