#  The Pollard-Rho algorithm for finding discrete logarithm in E(GF(p))

*(you can checkout branch `pollard-rho-Fp` for the algorithm in `GF(p)` instead)*

## Introduction

I'm back, with a **part III** of a trilogy that no one cares *(but I still did it anyways 🤩)*

This is an algorithm that tries to find the solution to the following problem: Given point `G` and `k*G` on the curve `y^2 = x^3 + ax + b` in `GF(p)` with a known multiplicative **PRIME** order `n`, find `k`.

This repo is made for educational purposes, striving for optimal time and space complexity in the discrete log for elliptic curves *(in short Weierstrass form)* in `GF(p)`.

Utilizes **GMP version 6.2.1**.

## Usages

```
git clone https://github.com/Mistsuu/BabyStepGiantStepC
cd BabyStepGiantStepC
git checkout pollard-rho
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
<n>
```

⚠️ It is required that `n` must be **PRIME** in order for the algorithm to work, else the algorithm just spits out `None`.

### Customize your number of threads

You can run `./dlog -t <num_threads>` to specify the number of threads used in multithreading part. If not specified, the default value for `num_threads` is `4`. 

### Limit your memory

You don't need to because in this version, the memory usage is very minimal *(less than `1MB` is expected for most applications)*

### Some custom parameters
- `-c <num_cache_items>` 
This value might allow you to have less errors processing data in a multithread setting. You can customize it to be any value `> 0` but it is advised to put it to any value from `4` to `10`. You can try change it higher to experiement with the runtime. However, the higher the value gets, the less the impact it would have on the runtime. Default value is `4`.

- `-r <nrandpoints>`
Set up the number of random points on the curve for the Pollard-Rho algorithm. You can experiment with this value, but any value from `20` and upper would work just as fine and won't affect much runtime. Default value is `20`.

### Output

The output will either be a number, or `None`, or some error data *(only happens in the case of memory error or thread creation error, which is not often as long as `n` is small enough that its square root fits 64-bits)*.

The program outputs `None` in the following cases:

- The program detects that no solution can be found.
- `n` is not a prime number or it is not positive.
- `num_cache_items` is `0`.
- `nrandpoints` is `< 2`. 

## Example

For example, to recover `k` from:
```
G = (437471552757133390 : 354281835126765881 : 1)
k*G = (295738136557598210 : 89525692852745488 : 1)
```
in curve `y^2 = x^3 + 1986076773346522069111732327339*x + 808177731529494834911895879646 in GF(13276420418771432419898581447951)` where we know the order of `G` is `857765763956341` using `8` threads.

We can run `./dlog -t 8` & supply the following input:
```
1986076773346522069111732327339
808177731529494834911895879646
13276420418771432419898581447951
12752653901711390718579996242468
9102988295173351464328400869432
6229151533029573525925181224628
1280290834308035816922484971565
857765763956341
```

Which gives the output:
```
690204827669615
```

If compiled with `BUILD=verbose` *(see the next section, **Compile modes**, for more detail)*, it will produce some outputs like this:
```
[debug] curve: 
[debug]    Elliptic Curve y^2 = x^3 + 1986076773346522069111732327339*x + 808177731529494834911895879646 in GF(13276420418771432419898581447951)
[debug] G: 
[debug]    (12752653901711390718579996242468 : 9102988295173351464328400869432 : 1)
[debug] kG: 
[debug]    (6229151533029573525925181224628 : 1280290834308035816922484971565 : 1)
[debug] G_mult_order = 857765763956341
[debug] n_threads = 4
[debug] n_cache_items = 10
[debug] n_rand_items = 20
[debug] index_size_limbs = 1
[debug] item_size_limbs = 2
[debug] Collision found!
[debug] Finished. Found k = 690204827669615
690204827669615
```

You can see some input examples provided in the `examples/` folder.


## Compile modes

Running `make`, you can specify `BUILD` variable to `release`, `verbose`, `memcheck`, `static` which creates different kind of builds:

- `release`: Using `dlog` will produce no debug output. *(default)*
- `verbose`: Using `dlog` will produce debug output.
- `memcheck`: Which just compiles the code with `-fsanitize=address`. Helpful in looking for memory leaks in the code.
- `static`: Creates a static version of `release` build of `dlog`.
- `allwarn`: If you want to compile with a lot of error displaying 🥰

You can also run `make libdlogefp` to build `so/libdlogefp.so` that you can use with the `so/dlog_EFp.py` *(but it gets a very weird bug that if you run the function `discrete_log_EFp()` a-lot, suddenly there's a bottleneck that cause the code to run much slower...)*

## Comparisons with the [parent project](https://github.com/Mistsuu/BabyStepGiantStepC)

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

  in the parent project from **34 seconds** to about **6 seconds** *(best running time, because this algorithm is very random, so sometimes it can get to **35 seconds**...)* using *4 threads* on **Intel(R) Core(TM) i5-10300H CPU @ 2.50GHz**. But on average, this should be faster than the **Baby Step Giant Step** one.

- Uses `O(1)` memory with a modest constant. *(so I could delete the memory handling code in satisfaction urg feeling so good)*

- The slowest operation's time complexity is `O(sqrt(n))`, where `n` is the order of `G`.

### Bads ❌ 

- `malloc()` fails are STILL handled by quick-and-dirty-`exit(-1)`s, (not sure if it's bad or not?)
- Runtime is not stable. The fastest ones are `4x` faster than the slowest ones due to the fact that the algorithm is probabilistic rather than deterministic like **Baby Step, Giant Step**.
- If my logic to detect that `kG` is not `k*G` using Weil Pairing is not enough, we might go to an infinite loop.

## How it works

I just give a link here because I know people can do it much better than me. The link I provided has very good formulas, explanations that can satisfy your satisfaction. The keyword to Google search is `Pollard-rho for discrete logarithm`.

- [Random Walks Revisited: Extensions of
  Pollard’s Rho Algorithm for Computing
  Multiple Discrete Logarithms](https://ac.informatik.uni-freiburg.de/publications/publications/sac01.pdf)

I use Teske's random function and apply Brent cycle-finding algorithm mentioned in this paper.

## Optimization

### Point addition optimization

I decided to represent points in the projective coordinate `(X:Y:Z)` so that no inversion mod `p` is required during computing new points *(they're really computationally expensive)*. I use these two formulas to add and double points in the curve:

- `add`: [addition-madd-1998-cmo](https://hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#addition-madd-1998-cmo) *(`Z2=1` because the algorithm only requires us to add with some fixed, random points)*
- `dbl`: [doubling-dbl-2007-bl](https://hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#doubling-dbl-2007-bl)

However, I regret it instantly because I still have to do the inverse thingy. The reason is that the algorithm requires you to have a function that uniquely maps a point to a certain value. 

### Reduce mod `p` optimization

I decided to put numbers in [Montgomery form](https://www.wikiwand.com/en/Montgomery_modular_multiplication) to do arithmetics with them. It helps me in a lot of ways.

Normally when you do `a*b mod p`, you have to do:

- `X = a*b`
- then `Y = X mod p`

However, the `mod` operation is so expensive that you can basically replace it with `3` multiplications and `1` subtraction and it will have less runtime. The Montgomery form allows you to do a map:

```
(a*R mod p), (b*R mod p) -> (a*b*R mod p)
```

where `R` is some random number you choose. While this map still requires you to do `mod`, but now it's in `mod R` instead of `mod p`. If you choose `R` to be `2^n` then `mod R` is just an `and` operation and that's how you save time baby 🤑🤑🤑!!! Better, if you choose `R` to be `mp_bits_per_limb` times the number of `mp_limb_t`s of `p`, you can just omit the first limbs :happy:

That's the first reason. The second reason is that after every arithmetic operations (`*`, `+`, `-`), the result always has the same number of `mp_limb_t`s as the inputs, which means that I don't have to keep track the number of limbs to allocate the right amount of memory. And also I can just write some Python code to automatically generate the C code for that part *(yey automation)*

And, because `(X:Y:Z)` is the same as `(XR:YR:ZR)`, so we can apply the arithmetics again and again to the points without having to transform the results to put it into the correct form.

### Multithreading

If we have `t` threads running at the same time applying the `Pollard-Rho` algorithm independently, we will have a speedup of `sqrt(t)` times. If the threads can communicate with each other, we will have a speed up of `t` times. That's what I (try) to do in this code, multi-threading and making them communicate with each other. 

Every thread will have a situation of *"one write, many read"*. It is a situation where one thread writes a value to a shared memory and the other threads will try to read it, compare with their data to get the result. 

I don't use locks, out of fear that it might hinder the finding process. Instead I decided to spread the writes into many memory slots so that the reading thread doesn't read the data at the same place of the writing thread writes. That's what the `-c` option are for, it specifies the number of those slots. 

While this mean that we might misread some values, it provide enough speedup so I just roll with this route :) *(and, yes, because I'm too lazy to implement other ways)* *(it probably also explain why we get such drastic fluctuation in runtime though...)*
