# Baby Step Giant Step for Elliptic Curves in F_p Parallelized, but is in C.

*(you can checkout `pollard-rho` branch for the Pollard-Rho implementation of the discrete log algorithm, which does the same thing as this but is **MUCH, MUCH faster** and consume **WAY, WAY much less memory**.)*

*(also `pollard-rho-Fp` for discrete log in `GF(p)` instead, `baby-Fp` is now for educational purposes only.)*

## Introduction

A child project spawned from [baby-giant-Fp-parallel](https://github.com/Mistsuu/baby-giant-Fp-parallel), and it's is in **C**.

This is an algorithm that tries to find the solution to the following problem: Given point `G` and `k*G` on the curve `y^2 = x^3 + ax + b` in `GF(p)` and an upper bound `upper_k` of `k`, find `k`.

This repo is made for educational purposes, striving for optimal time and space complexity in baby step & giant step algorithm for elliptic curves in `GF(p)`.

Utilizes **GMP version 6.2.1**.

## Usages

```
git clone https://github.com/Mistsuu/dlogC
cd dlogC
make -j16
```
to produce `./dlog`.

### Input supplying

To use `./dlog`, we supply input as a list of numbers seperated by a newline (`\n`) in the following format:

```
<curve.a>
<curve.b>
<curve.p>
<X(G)>
<Y(G)>
<X(kG)>
<Y(kG)>
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
```
G = (437471552757133390 : 354281835126765881 : 1)
k*G = (295738136557598210 : 89525692852745488 : 1)
```
in curve `y^2 = x^3 + 448019786388741247*x + 544225411105375163 in GF(593010448435692503)` where we know the order of `G` is `593010448361862286` using `8` threads limited to **4GB** of memory.

We can run `./dlog 8 4G` & supply the following input:
```
448019786388741247
544225411105375163
593010448435692503
437471552757133390
354281835126765881
295738136557598210
89525692852745488
593010448361862286
```

Which gives the output:
```
234176126564864674
```

If compiled with `BUILD=verbose` *(see the next section, **Compile modes**, for more detail)*, it will produce some outputs like this:
```
[debug] curve: 
[debug]    Elliptic Curve y^2 = x^3 + 448019786388741247*x + 544225411105375163 in GF(593010448435692503)
[debug] G: 
[debug]    (437471552757133390 : 354281835126765881 : 1)
[debug] kG: 
[debug]    (295738136557598210 : 89525692852745488 : 1)
[debug] upper_k = 593010448361862286
[debug] memory limit: 4294967296 bytes = 4096.000000 MB = 4.000000 GB
[debug] n_threads = 8
[debug] index_size_bytes = 4
[debug] item_size_bytes = 8
[debug] index_size_limbs = 1
[debug] item_size_limbs = 1
[debug] n_partitions = 19
[debug] n_items = 178956969
[debug] size buffer: 4294967280 bytes = 4095.999985 MB = 4.000000 GB

[debug] === Running partition 0 === (found k value in this partition will be added with 32025596753666961*0 to get the actual k)
[debug] Filling L buffer...
[debug] Filling L took 20.307818 seconds.
[debug] Filling R buffer...
[debug] Filling R took 21.101883 seconds.
[debug] Sorting L buffer...
[debug] Sorting L took 62.114992 seconds.
[debug] Sorting R buffer...
[debug] Sorting R took 62.804649 seconds.
[debug] Searching in L & R buffers...
[debug] Searching took 4.806363 seconds.
[debug] Cannot search for equal values in L & R buffers!

[debug] === Running partition 1 === (found k value in this partition will be added with 32025596753666961*1 to get the actual k)
[debug] Filling R buffer...
[debug] Filling R took 20.896676 seconds.
[debug] Sorting R buffer...
[debug] Sorting R took 60.835373 seconds.
[debug] Searching in L & R buffers...
[debug] Searching took 4.849460 seconds.
[debug] Cannot search for equal values in L & R buffers!

[debug] === Running partition 2 === (found k value in this partition will be added with 32025596753666961*2 to get the actual k)
[debug] Filling R buffer...
[debug] Filling R took 22.212721 seconds.
[debug] Sorting R buffer...
[debug] Sorting R took 62.251165 seconds.
[debug] Searching in L & R buffers...
[debug] Searching took 4.489302 seconds.
[debug] Cannot search for equal values in L & R buffers!

[debug] === Running partition 3 === (found k value in this partition will be added with 32025596753666961*3 to get the actual k)
[debug] Filling R buffer...
[debug] Filling R took 20.863485 seconds.
[debug] Sorting R buffer...
[debug] Sorting R took 61.216875 seconds.
[debug] Searching in L & R buffers...
[debug] Searching took 4.551580 seconds.
[debug] Cannot search for equal values in L & R buffers!

[debug] === Running partition 4 === (found k value in this partition will be added with 32025596753666961*4 to get the actual k)
[debug] Filling R buffer...
[debug] Filling R took 20.377466 seconds.
[debug] Sorting R buffer...
[debug] Sorting R took 61.615065 seconds.
[debug] Searching in L & R buffers...
[debug] Searching took 4.581498 seconds.
[debug] Cannot search for equal values in L & R buffers!

[debug] === Running partition 5 === (found k value in this partition will be added with 32025596753666961*5 to get the actual k)
[debug] Filling R buffer...
[debug] Filling R took 20.770146 seconds.
[debug] Sorting R buffer...
[debug] Sorting R took 61.968580 seconds.
[debug] Searching in L & R buffers...
[debug] Searching took 4.491920 seconds.
[debug] Cannot search for equal values in L & R buffers!

[debug] === Running partition 6 === (found k value in this partition will be added with 32025596753666961*6 to get the actual k)
[debug] Filling R buffer...
[debug] Filling R took 20.651545 seconds.
[debug] Sorting R buffer...
[debug] Sorting R took 61.693792 seconds.
[debug] Searching in L & R buffers...
[debug] Searching took 4.452649 seconds.
[debug] Cannot search for equal values in L & R buffers!

[debug] === Running partition 7 === (found k value in this partition will be added with 32025596753666961*7 to get the actual k)
[debug] Filling R buffer...
[debug] Filling R took 20.949573 seconds.
[debug] Sorting R buffer...
[debug] Sorting R took 61.580475 seconds.
[debug] Searching in L & R buffers...
[debug] Searching took 1.156795 seconds.
[debug] Found k = 9996949289195947.
234176126564864674
```

You can see some input examples provided in the `examples/` folder.


## Compile modes

Running `make`, you can specify `BUILD` variable to `release`, `verbose`, `memcheck`, `static` which creates different kind of builds:

- `release`: Using `dlog` will produce no debug output. *(default)*
- `verbose`: Using `dlog` will produce debug output such as:
  - The size of allocated memory to construct the `L` and `R` `char` buffers in the baby step giant step algorithm.
  - Time took for each sub-operations.
  - And many more...
- `memcheck`: Which just compiles the code with `-fsanitize=address`. Helpful in looking for memory leaks in the code.
- `static`: Creates a static version of `release` build.

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

  in the parent project from **10 minutes** to **34 seconds** using *8 threads* on **Intel(R) Core(TM) i5-10300H CPU @ 2.50GHz**.

- Also reduces some more memory.

### Bads ❌ 

- Cannot parallelize array-sorting in the 2nd step without causing some unknown bottleneck problem.
- `malloc()` fails are handled by quick-and-dirty-`exit(-1)`s, (not sure if it's bad or not?)
- Code is probably unnecessarily long because of some function is built such as `xmul()`s that is not used at all. It's only there because I just like to implement them for no other reasons than curiosity.

## How it works

### Baby Step Giant Step Algorithm Basic

The algorithm does this by storing `n+1` points (`n = isqrt(upper_k)`) to 2 arrays: `L` and `R`:

- `L` stores `0*G`, `1*G`, `2*G`, ..., `n*G`.
- `R` stores `k*G`, `(k-n)*G`, ..., `(k-n*n)*G`.

If we can find `l*G` in `L` and `(k-r*n)*G` in `R` that `l*G == (k-r*n)*G`, we can solve for `k = l + r*n`.

### Sub-operations

This code divides the process into 3 sub-operations:

#### Fill `L` & `R` 

Filling `L` and `R` with the above points. This part can be space-optimized by storing each point like this:

```
+---------+-------------------------+
|  index  |        X(index*G)       |
+---------+-------------------------+
```

We store `X` coordinates of the points only. Because if `X(l*G) == X((k-r*n)*G)`, we can still recover `k` from `k = l + r*n`, or `k = r*n - l`.

I use multi-threading in this sub-operations to speed up the filling.

Also I choose to represent point's coordinates in `X/Z` form, so that it only takes one inversion for each point to speed up this part even more :>

#### Sort `L` & `R`

Sort `L` and `R` by each element's `X` coordinate so that we can search for equal values in `O(N)` time. *(this operation takes `O(NlogN)`, however)*, where `N` is `isqrt(upper_k)`.

It uses **Quick Sort** to sort the array, allowing an in-place memory sort, thus requiring no additional memory usages.

#### Search `L` and `R`

After the arrays are sorted, we search the equal `X` values in them, then deduce `k` to get the result.
