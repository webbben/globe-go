# Globe-Go: ASCII globe generator, ported from Rust to Go

This is a Go wrapper around a Rust ASCII globe generator. Originally written by adamsky, here I've made it so the main pieces of the globe functionality can be used in Go applications.

## How does it work?

I've added a few FFI-compatible (Foreign Function Interface) external functions to the rust code, so we can build the rust project into a binary, and then in Go we can use `cgo` to access these functions.

From the Go side, it's the same process as calling C code.

## Build

Follow these steps to build this source code on your own machine.

### Prerequisites

To build the Rust code, you need to have the [Rust programming language installed](https://www.rust-lang.org/tools/install).

To build the Go code, similarly, you will need Go installed. But I'm assuming anyone who wants to use this Go wrapper in the first place already has this covered.

#### Cross-Platform Builds

You may need to install different things to be able to build for other platforms, such as if you are on a mac but want to build for linux.

First, you need to install the right target:

```shell
# build target for linux
rustup target add x86_64-unknown-linux-gnu

# build target for windows
rustup target add x86_64-pc-windows-gnu
```

The rust compiler will probably tell you if you are missing this, and instruct you to run one of these commands.

#### ARM64 / M1 Architecture

Apparently compiling to other platforms from a macbook with M1 architecture is pretty difficult, and I haven't succeeded so far. I challenge someone to figure out the prerequisites and/or build script and add information here.

### Build Steps

1. First, build the Rust code. Since each OS will have a different type of binary it wants, choose the `build_<os>.sh` file matching your OS.

```shell
cd <project-root>
bash build_darwin.sh
```

The build artifact will be in `/bin` (from the project root directory). But you won't need to manually access this - the Go code already knows to check here for the OS-specific binary.

2. Build the Go module. You can just do `go build` from the project root directory.

## Usage

> TODO
