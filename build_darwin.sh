#!/bin/bash

cd rust-src
cargo build -p globe --release
cp target/release/libglobe.dylib ../bin/libglobe.dylib