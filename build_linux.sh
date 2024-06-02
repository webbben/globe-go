#!/bin/bash

# make sure the target is installed
rustup target add x86_64-unknown-linux-gnu

cd rust-src
cargo build -p globe --release --target x86_64-unknown-linux-gnu
cp target/x86_64-unknown-linux-gnu/release/libglobe.so ../bin/libglobe.so