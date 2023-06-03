#ifndef EXMPN_H
#define EXMPN_H

#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

void       mpn_printf   (mp_limb_t* in, mp_size_t len);
mp_limb_t* mpn_init_cpyz(mpz_t x, size_t n);
mp_limb_t* mpn_init_zero(size_t n);
void       mpn2bytes    (unsigned char *str, mp_size_t len, const mp_limb_t *s1p, mp_size_t s1n);

void mpn_montgomery_mulmod_n(
    mp_limb_t* rp, 
    const mp_limb_t* s1p, const mp_limb_t* s2p, 
    const mp_limb_t* dp, const mp_limb_t* Dp, 
    mp_size_t n, 
    mp_limb_t* tp
);

void mpn_montgomery_sqrmod_n(
    mp_limb_t* rp, 
    const mp_limb_t* s1p, 
    const mp_limb_t* dp, const mp_limb_t* Dp,
    mp_size_t n, 
    mp_limb_t* tp
);

#define mpn_montgomery_addmod_n(rp, s1p, s2p, dp, n)            \
do {                                                            \
    if (mpn_add_n(rp, s1p, s2p, n) || mpn_cmp(rp, dp, n) >= 0)  \
        mpn_sub_n(rp, rp, dp, n);                               \
} while(0)

#define mpn_montgomery_lshift1mod_n(rp, s1p, dp, n)                \
do {                                                               \
    if (mpn_lshift(rp, s1p, n, 1) || mpn_cmp(rp, dp, n) >= 0)      \
        mpn_sub_n(rp, rp, dp, n);                                  \
} while(0)

#define mpn_montgomery_submod_n(rp, s1p, s2p, dp, n) \
do {                                                 \
    if (mpn_sub_n(rp, s1p, s2p, n))                  \
        mpn_add_n(rp, rp, dp, n);                    \
} while(0)

#endif