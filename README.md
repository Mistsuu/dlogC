# Baby Step Giant Step for Elliptic Curves in F_p Parallelized, but is in C.

*(you can checkout branch `baby-Fp` for baby-step-giant-step algorithm in `Fp` instead)*

## Introduction

A child project of a child project spawned from [baby-giant-Fp-parallel](https://github.com/Mistsuu/baby-giant-Fp-parallel), and it's is in **C**.

This is an algorithm that tries to find the solution to the following problem: Given `G` and `G^k mod p`, find `k`.

This repo is made for educational purposes, and for fun.

Utilizes **GMP version 6.2.1**.

## Usages

```
git clone https://github.com/Mistsuu/BabyStepGiantStepC
cd BabyStepGiantStepC
make -j16
```

to produce `./dlog`.

### Input supplying

To use `./dlog`, we supply input as a list of numbers seperated by a newline (`\n`) in the following format:

```
<p>
<G>
<G^k mod p>
<upper k bound>
```

### Customize your number of threads

You can run `./dlog <num_threads>` to specify the number of threads used in multithreading part. If not specified, the default value for `num_threads` is `4`. 

### Limit your memory

In the newer version, running `./dlog <num_threads> <memory_limit>` will run `./dlog` with limited memory. Suitable for machines where not enough RAM is available. The format of `<memory_limit>` is a `double` value followed by `M` - **megabytes** (*for example,* `4096M`) or `G` - **gigabytes** (*for example,* `1.25G`). If only a `double` value is given, the unit will be `G` by default.

If you don't specify this part at all, the program will use any amount of memory as it likes.

### Output

The output will either be a number *(a negative one is normal)*, or `None`, or some error data *(only happens in the case of memory error or thread creation error, which is not often as long as `<upper k bound>` is small enough that its square root fits 64-bits)*.

## Example

For example, to recover `k` from:

```python
G = 14484073937158643
pow(G, k, p) = 14484073937158643
p = 17524524814163177
```

where we know the order of `G` is `17524524814163176` using `8` threads limited to **1GB** of memory.

We can run `./dlog 8 1G` & supply the following input:

```
17524524814163177
14484073937158643
16341616837220933
17524524814163176
```

Which gives the output:

```
8766935282215935
```

If compiled with `BUILD=verbose` *(see the next section, **Compile modes**, for more detail)*, it will produce some outputs like this:

```
[debug] p: 
[debug]    17524524814163177
[debug] G: 
[debug]    14484073937158643
[debug] kG: 
[debug]    16341616837220933
[debug] upper_k = 17524524814163176
[debug] memory limit: 1073741824 bytes = 1024.000000 MB = 1.000000 GB
[debug] n_threads = 8
[debug] index_size_bytes = 4
[debug] item_size_bytes = 7
[debug] index_size_limbs = 1
[debug] item_size_limbs = 1
[debug] n_partitions = 8
[debug] n_items / partition = 48806445
[debug] size buffer: 1073741812 bytes = 1023.999989 MB = 1.000000 GB

[debug] === Running partition 1/8 === (found k value in this partition will be added with 2382069073538025*0 to get the actual k)
[debug] Filling L buffer...
[debug] Filling L took 0.690344 seconds.
[debug] Filling R buffer...
[debug] Filling R took 0.748650 seconds.
[debug] Sorting L buffer...
[debug] Sorting L took 15.774451 seconds.
[debug] Sorting R buffer...
[debug] Sorting R took 15.809021 seconds.
[debug] Searching in L & R buffers...
[debug] Searching took 1.191961 seconds.
[debug] Cannot search for equal values in L & R buffers!

[debug] === Running partition 2/8 === (found k value in this partition will be added with 2382069073538025*1 to get the actual k)
[debug] Filling R buffer...
[debug] Filling R took 0.771460 seconds.
[debug] Sorting R buffer...
[debug] Sorting R took 16.296818 seconds.
[debug] Searching in L & R buffers...
[debug] Searching took 1.207998 seconds.
[debug] Cannot search for equal values in L & R buffers!

[debug] === Running partition 3/8 === (found k value in this partition will be added with 2382069073538025*2 to get the actual k)
[debug] Filling R buffer...
[debug] Filling R took 0.725180 seconds.
[debug] Sorting R buffer...
[debug] Sorting R took 15.734189 seconds.
[debug] Searching in L & R buffers...
[debug] Searching took 1.193694 seconds.
[debug] Cannot search for equal values in L & R buffers!

[debug] === Running partition 4/8 === (found k value in this partition will be added with 2382069073538025*3 to get the actual k)
[debug] Filling R buffer...
[debug] Filling R took 0.708008 seconds.
[debug] Sorting R buffer...
[debug] Sorting R took 16.147236 seconds.
[debug] Searching in L & R buffers...
[debug] Searching took 1.195337 seconds.
[debug] Found k = 1620728061601860.
8766935282215935
```

You can see some input examples provided in the `examples/` folder.


## Compile modes

Running `make`, you can specify `BUILD` variable to `release`, `verbose`, `memcheck` which creates different kind of builds:

- `release`: Using `dlog` will produce no debug output. *(default)*
- `verbose`: Using `dlog` will produce debug output such as:
  - The size of allocated memory to construct the `L` and `R` `char` buffers in the baby step giant step algorithm.
  - Time took for each sub-operations.
  - And many more...
- `memcheck`: Which just compiles the code with `-fsanitize=address`. Helpful in looking for memory leaks in the code.
- `static`: Provides static build for the code.

## How it works

### Baby Step Giant Step Algorithm Basic

The algorithm does this by storing `n+1` points (`n = isqrt(G.order())`) to 2 arrays: `L` and `R`:

- `L` stores `G^0`, `G^1`, `G^2`, ..., `G^n`.
- `R` stores `G^k`, `G^(k-n)`, ..., `G^(k-n*n)`.

If we can find `G^l` in `L` and `G^(k-r*n)` in `R` that `G^l == G^(k-r*n)`, we can solve for `k = l + r*n`.

### Sub-operations

This code divides the process into 3 sub-operations:

#### Fill `L` & `R` 

Filling `L` and `R` with the mentioned values.

```
+---------+-------------------------+
|  index  |          G^index        |
+---------+-------------------------+
```

I use multi-threading in this sub-operations to speed up the filling.

#### Sort `L` & `R`

Sort `L` and `R` by each element's `X` coordinate so that we can search for equal values in `O(N)` time. *(this operation takes `O(NlogN)`, however)*

It uses **Quick Sort** to sort the array, allowing an in-place memory sort, thus requiring no additional memory usages.

#### Search `L` and `R`

After the arrays are sorted, we search the equal `X` values in them, then deduce `k` to get the result.