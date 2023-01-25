# BabyStepGiantStepC
Same thing as BabyStepGiantStep, but is in C.
This repo is for educational purposes, striving for optimal time and space complexity in baby step & giant step algorithm for elliptic curve's.

Compiled and utilized using GMP version 6.2.1.

## Goods ✅
- Reduces runtime for example test-case in Python from **10 minutes** to **37 seconds** using *4-cores*.
- Also reduces some more memory compared to Python code.

## Bads ❌ 
- Cannot parallelize array-sorting in the 2nd step without causing some unknown bottleneck problem.
- `malloc()` fails are handled by quick-and-dirty-`exit(-1)`s, (not sure if it's bad or not?)
- Code is probably unecessarily long because of `xmul()`s or `xadd()`s and `xdbl()`s that are not used at all. To be honest, I just like to implement them while doing this, so that's the reason they are there...

## Running
```
git clone https://github.com/Mistsuu/BabyStepGiantStepC
cd BabyStepGiantStepC
make main
```
to compile `src/main.c` to `./main` in `BabyStepGiantStepC` folder.

Modify `src/main.c` following by `make` for a quick use.