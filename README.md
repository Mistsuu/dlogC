# The Pollard-Rho algorithm for finding discrete logarithm in GF(p)

*(if no solution exists, this program wil run foreverrrrrrr)*

## Introduction

Part **III.75** of the [trilogy](https://github.com/Mistsuu/BabyStepGiantStepC).

This is an algorithm that tries to find the solution to the following problem: Given `G` and `G^k mod p` with a known multiplicative **PRIME** order `n`, find `k`.

This repo is made for educational purposes, and for fun.

Utilizes **GMP version 6.2.1**.

## Usages

```
git clone https://github.com/Mistsuu/BabyStepGiantStepC
cd BabyStepGiantStepC
git checkout pollard-rho-Fp
make -j16
```

to produce `./dlog`.

### Input supplying

To use `./dlog`, we supply input as a list of numbers seperated by a newline (`\n`) in the following format:

```
<p>
<G>
<G^k mod p>
<n>
```

⚠️ It is required that `n` must be **PRIME** in order for the algorithm to work, else the algorithm just spits out `None`.

### Customize your number of threads

You can run `./dlog -t <num_threads>` to specify the number of threads used in multithreading part. If not specified, the default value for `num_threads` is `4`. 

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

- The program detects that no solution can be found *(but in a very very small case)*
- `n` is not a prime number or it is not positive.
- `num_cache_items` is `0`.
- `nrandpoints` is `< 2`. 

## Example

For example, to recover `k` from:

```python
G = 63848476900962761921955720089397765658257438572112067334361113887001094555500037
pow(G, k, p) = 124598459814481000866971099140174848414430762122941331126292491764551683720313014
p = 17524524814163177
```

where we know the order of the multiplicative group generated by `G` is `2259283924057529` using `4` threads.

We can run `./dlog -t 4` & supply the following input:

```
140601987068392810555209463610249011351761931006927193734529493874226799329006719
63848476900962761921955720089397765658257438572112067334361113887001094555500037
124598459814481000866971099140174848414430762122941331126292491764551683720313014
2259283924057529
```

Which gives the output:

```
814905979740757
```

If compiled with `BUILD=verbose` *(see the next section, **Compile modes**, for more detail)*, it will produce some outputs like this:

```
[debug] p: 
[debug]    140601987068392810555209463610249011351761931006927193734529493874226799329006719
[debug] G: 
[debug]    63848476900962761921955720089397765658257438572112067334361113887001094555500037
[debug] kG: 
[debug]    124598459814481000866971099140174848414430762122941331126292491764551683720313014
[debug] G_mult_order = 2259283924057529
[debug] n_threads = 4
[debug] n_cache_items = 10
[debug] n_rand_items = 20
[debug] index_size_limbs = 1
[debug] item_size_limbs = 5
[debug] Collision found!
[debug] Finished. Found k = 814905979740757
814905979740757
```

You can see some input examples provided in the `examples/` folder.


## Compile modes

Running `make`, you can specify `BUILD` variable to `release`, `verbose`, `memcheck`, `static` which creates different kind of builds:

- `release`: Using `dlog` will produce no debug output. *(default)*
- `verbose`: Using `dlog` will produce debug output.
- `memcheck`: Which just compiles the code with `-fsanitize=address`. Helpful in looking for memory leaks in the code.
- `static`: Creates a static version of `release` build of `dlog`.
- `allwarn`: If you want to compile with a lot of error displaying 🥰

You can also run `make libdlogfp` to build `so/libdlogfp.so` that you can use with the `so/dlog_Fp.py` *(but it gets a very weird bug that if you run the function `discrete_log_Fp()` a-lot, suddenly there's a bottleneck that cause the code to run much slower...)*

## How it works

I just give a link here because I know people can do it much better than me. The link I provided has very good formulas, explanations that can satisfy your satisfaction. The keyword to Google search is `Pollard-rho for discrete logarithm`.

- [Random Walks Revisited: Extensions of
  Pollard’s Rho Algorithm for Computing
  Multiple Discrete Logarithms](https://ac.informatik.uni-freiburg.de/publications/publications/sac01.pdf)

I use Teske's random function and apply Brent cycle-finding algorithm mentioned in this paper.

## Optimization

*(i'm lazy so I just copy the description from the `pollard-rho` branch... please don't hit me if i'm wrong some part :()*

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

### Multithreading

If we have `t` threads running at the same time applying the `Pollard-Rho` algorithm independently, we will have a speedup of `sqrt(t)` times. If the threads can communicate with each other, we will have a speed up of `t` times. That's what I (try) to do in this code, multi-threading and making them communicate with each other. 

Every thread will have a situation of *"one write, many read"*. It is a situation where one thread writes a value to a shared memory and the other threads will try to read it, compare with their data to get the result. 

I don't use locks, out of fear that it might hinder the finding process. Instead I decided to spread the writes into many memory slots so that the reading thread doesn't read the data at the same place of the writing thread writes. That's what the `-c` option are for, it specifies the number of those slots. 

While this mean that we might misread some values, it provide enough speedup so I just roll with this route :) *(and, yes, because I'm too lazy to implement other ways)* *(it probably also explain why we get such drastic fluctuation in runtime though...)*
