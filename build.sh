#!/usr/bin/env bash
# Builds the wasm64 core. Needs nightly: wasm64-unknown-unknown is a Tier 3
# target with no precompiled core, and core::arch::wasm64 is still unstable.
#
# Note the absence of +atomics / --shared-memory: the memory is NOT shared, so
# the page needs no COOP/COEP headers and works on plain static hosting.
set -euo pipefail
rustup toolchain install nightly --profile minimal --component rust-src

RUSTFLAGS="-C target-feature=+bulk-memory,+mutable-globals \
 -C link-arg=--import-memory \
 -C link-arg=--max-memory=17179869184 \
 -C link-arg=--stack-first -C link-arg=-zstack-size=1048576 \
 -C link-arg=--export=__heap_base -C link-arg=--export=__data_end" \
cargo +nightly build --release --target wasm64-unknown-unknown \
  -Z build-std=core,compiler_builtins

cp target/wasm64-unknown-unknown/release/fqclean.wasm ./fqclean.wasm
echo "built: $(stat -c%s fqclean.wasm) bytes"
