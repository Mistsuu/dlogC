#include "ex_mpn.h"
#include "mem.h"

// =================================================================================
//                                 MEMORY STUFFS
// =================================================================================

void mpn_printf(mp_limb_t* mpn_in, mp_size_t mpn_len)
{
    mpz_t mpz_in;
    mpz_in->_mp_alloc = mpn_len;
    mpz_in->_mp_size  = mpn_len;
    mpz_in->_mp_d     = mpn_in;
    mpz_out_str(stdout, 10, mpz_in);
}

mp_limb_t* mpn_init_cpyz(mpz_t x, size_t n)
{
    const mp_limb_t* x_limbs = mpz_limbs_read(x);
    size_t           x_size  = mpz_size(x);

    mp_limb_t* x_copied_limbs = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
    memset(x_copied_limbs, 0, sizeof(mp_limb_t) * n);
    memcpy(x_copied_limbs, x_limbs, sizeof(mp_limb_t) * (x_size > n ? n : x_size));
    return x_copied_limbs;
}

mp_limb_t* mpn_init_zero(size_t n)
{
    mp_limb_t* x_limbs = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
    memset(x_limbs, 0, sizeof(mp_limb_t) * n);
    return x_limbs;
}

void mpn_cpyz(mp_limb_t* rop, mpz_t x, size_t n)
{
    const mp_limb_t* x_limbs = mpz_limbs_read(x);
    size_t           x_size  = mpz_size(x);
    memset(rop, 0, sizeof(mp_limb_t) * n);
    memcpy(rop, x_limbs, sizeof(mp_limb_t) * (x_size > n ? n : x_size));
}

void mpn2bytes(unsigned char *str, mp_size_t len, const mp_limb_t *s1p, mp_size_t s1n)
{
    // If I don't do this, 
    // mpn_get_str() will return some crazy stuffs :'<
    if (!s1n || mpn_zero_p(s1p, s1n)) {
        memset(str, 0, len);
        return;
    }
    
    // According to the doc,
    // "The most significant limb of the input {s1p, s1n} must be non-zero."
    // Ignore it and you'll get a very nasty SIGSEGV :(
    mp_size_t actual_len = 0;
    while (s1n > 0 && !s1p[s1n - 1])
        s1n--;
    if (s1n > 0)
        actual_len = mpn_get_str(str, 256, (mp_limb_t *)s1p, s1n);

    // todo: i'm glad that i didn't delete this shit
    if (actual_len > len) {
        printf("[error] omg plz, whyyyyy mpn_get_str cannot fit in mpn2bytes: actual_len=%ld, len=%ld???????\n", actual_len, len);
        printf("[error] s1p: \n");
        for (mp_size_t i = 0; i < s1n; ++i) {
            printf("       0x%016lx\n", s1p[i]);
        }
        exit(-1);
    }
    
    memmove(str+len-actual_len, str, actual_len);
    memset(str, 0, len-actual_len);
}

// =================================================================================
//                                 ARITHMETICS STUFFS
// =================================================================================

/*
    mpn_montgomery_mulmod_n:
        This function performs a map:
            (s1 * R mod d) * (s2 * R mod d) -> (s1*s2 * R mod d)
        where R is a number satisfying:
            - is 2^(n * mp_bits_per_limb) 
            - gcd(R, d) == 1 (just put an odd number)

        Argument size should be:
            - {s1p, n}
            - {s2p, n}
            - {dp, n}
            - {Dp, n} = -{dp, n}^-1 mod R
            - {rp, n} (return value)
            - {tp, 6*n} (for temporary placement)
*/
void mpn_montgomery_mulmod_n(
    mp_limb_t* rp, 
    const mp_limb_t* s1p, const mp_limb_t* s2p, 
    const mp_limb_t* dp, const mp_limb_t* Dp, 
    mp_size_t n, 
    mp_limb_t* tp
)
{
    mp_limb_t* t1p = tp;
    mp_limb_t* t2p = &t1p[2*n];
    mp_limb_t* t3p = &t2p[2*n];
    mpn_mul_n(t1p, s1p, s2p, n);
    mpn_mul_n(t2p, t1p, Dp,  n);
    mpn_mul_n(t3p, t2p, dp,  n);
    if (mpn_add_n(t3p, t3p, t1p, 2*n) || mpn_cmp(&t3p[n], dp, n) >= 0)
        mpn_sub_n(rp, &t3p[n], dp, n);
    else
        mpn_copyd(rp, &t3p[n], n);
}

/*
    mpn_montgomery_sqrmod_n:
        This function performs a map:
            (s1 * R mod d) * (s1 * R mod d) -> (s1*s1 * R mod d)
        where R is a number satisfying:
            - is 2^(n * mp_bits_per_limb) 
            - gcd(R, d) == 1 (just put an odd number)

        Argument size should be:
            - {s1p, n}
            - {dp, n}
            - {Dp, n} = -{dp, n}^-1 mod R
            - {rp, n} (return value)
            - {tp, 6*n} (for temporary placement)
*/
void mpn_montgomery_sqrmod_n(
    mp_limb_t* rp, 
    const mp_limb_t* s1p, 
    const mp_limb_t* dp, const mp_limb_t* Dp,
    mp_size_t n, 
    mp_limb_t* tp
)
{
    mp_limb_t* t1p = tp;
    mp_limb_t* t2p = &t1p[2*n];
    mp_limb_t* t3p = &t2p[2*n];
    mpn_sqr  (t1p, s1p, n);
    mpn_mul_n(t2p, t1p, Dp,  n);
    mpn_mul_n(t3p, t2p, dp,  n);
    if (mpn_add_n(t3p, t3p, t1p, 2*n) || mpn_cmp(&t3p[n], dp, n) >= 0)
        mpn_sub_n(rp, &t3p[n], dp, n);
    else
        mpn_copyd(rp, &t3p[n], n);
}