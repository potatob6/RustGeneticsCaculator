## Rust Genetics Caculator
---
### Usage
&emsp;&emsp;A plants gene caculator for "Rust" game. The program use greedy algorithm to approach local optimal crossbreed solution.

&emsp;&emsp;The program will output the crossbreed steps from exists gene to target gene. (e.g. #1 and #2 is serial step. #2-1 and #2-2 the parallel steps that means can crossbreed at the same time and they will not depent on both).

&emsp;&emsp;The program only support Chinese now, maybe it will support English later.

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