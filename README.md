# Baby Step Giant Step for Elliptic Curves in F_p Parallelized, but is in C.

## Introduction

A child project spawned from [baby-giant-Fp-parallel](https://github.com/Mistsuu/baby-giant-Fp-parallel), and it's is in **C**.

This is an algorithm that tries to find the solution to the following problem: Given point `G` and `k*G` on the curve `y^2 = x^3 + ax + b` in `GF(p)`, find `k`.

This repo is made for educational purposes, striving for optimal time and space complexity in baby step & giant step algorithm for elliptic curves in `GF(p)`.

Utilizes **GMP version 6.2.1**.

## Running

```
git clone https://github.com/Mistsuu/BabyStepGiantStepC
cd BabyStepGiantStepC
make
```

to compile `src/main.c` to `./main` in `BabyStepGiantStepC` folder.

Modify `src/main.c` following by a `make` for a quick use.

## Makefile compile modes

In `Makefile`, you can modify `BUILD` variable to `release`, `verbose`, `memcheck` which creates different kind of builds:

- `release`: Using `dlog` will produce no debug output.
- `verbose`: Using `dlog` will produce debug output such as:
  - The size of allocated memory to construct the `L` and `R` `char` buffers in the baby step giant step algorithm.
  - Time took for each sub-operations.
- `memcheck`: Which just compiles the code with `-fsanitize=address`, that checks for memory leaks, but it's got a weird bug that saying I freed the wrong stuffs, so idk...

## Comparisons with the [parent project](https://github.com/Mistsuu/baby-giant-Fp-parallel)

### Goods ✅

- Reduces runtime for the **example test-case**:

  ```
  curve: y^2 = x^3 + ax + b in GF(p)
  where:
  	a = 1986076773346522069111732327339
  	b = 808177731529494834911895879646
  	p = 13276420418771432419898581447951
  	
  generator_point = (12752653901711390718579996242468:9102988295173351464328400869432:1)
  generator_order = 857765763956341
  ```

  in the parent project from **10 minutes** to **35 seconds** using *4-cores* on **Intel(R) Core(TM) i5-10300H CPU @ 2.50GHz**.

- Also reduces some more memory.

### Bads ❌ 

- Cannot parallelize array-sorting in the 2nd step without causing some unknown bottleneck problem.
- `malloc()` fails are handled by quick-and-dirty-`exit(-1)`s, (not sure if it's bad or not?)
- Code is probably unnecessarily long because of `xmul()`s or `xadd()`s and `xdbl()`s that are not used at all. To be honest, I just like to implement them while doing this, so that's the reason they are there...

## How it works

### Baby Step Giant Step Algorithm Basic

The algorithm does this by storing `n+1` points (`n = isqrt(G.order())`) to 2 arrays: `L` and `R`:

- `L` stores `0*G`, `1*G`, `2*G`, ..., `n*G`.
- `R` stores `k*G`, `(k-n)*G`, ..., `(k-n*n)*G`.

If we can find `l*G` in `L` and `(k-r*n)*G` in `R` that `l*G == (k-r*n)*G`, we can solve for `k = l + r*n`.

### Sub-operations

This code divides the process into 3 sub-operations:

#### Fill `L` & `R` 

Filling `L` and `R` with the above points. This part can be space-optimized by storing each point like this:

```
+---------+-------------------------+
|  index  |        X(index*G)		|
+---------+-------------------------+
```

We store `X` coordinates of the points only. Because if `X(l*G) == X((k-r*n)*G)`, we can still recover `k` from `k = l + r*n`, or `k = r*n - l `.

I use multi-threading in this sub-operations to speed up the filling.

Also I choose to represent point's coordinates in `X/Z` form, so that it only takes one inversion for each point to speed up this part even more :>

#### Sort `L` & `R`

Sort `L` and `R` by each element's `X` coordinate so that we can search for equal values in `O(N)` time. *(this operation takes `O(NlogN)`, however)*

It uses **Quick Sort** to sort the array, allowing an in-place memory sort, thus requiring no additional memory usages.

#### Search `L` and `R`

After the arrays are sorted, we search the equal `X` values in them, then deduce `k` to get the result.