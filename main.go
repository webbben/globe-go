package main

/*
#cgo LDFLAGS: -L./rust-src/target/release -lglobe
#include <stdio.h>
#include <stdlib.h>

// Declare the Rust function here
void greet();
*/
import "C"
import "fmt"

func main() {
	fmt.Println("Calling Rust function from Go:")
	C.greet()
}
