## Rust Genetics Caculator
---
### Usage
&emsp;&emsp;A plant's gene calculator for "Rust" game. The program uses greedy algorithm to approach local optimal crossbreed solution.

&emsp;&emsp;The program will output the crossbreed steps from existing gene to target gene. (e.g. #1 and #2 are serial steps. #2-1 and #2-2 are the parallel steps that means can crossbreed at the same time and they will not depend on both).

&emsp;&emsp;The program only supports Chinese now, maybe it will support English later.

---
### Build 
Rust project dependences:
- bigdecimal = "0.4.6"
- nohash = "0.2.0"

Command to build:
``` shell
cargo build -r
cargo run -r
```