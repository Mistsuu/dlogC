#ifndef NUM_H
#define NUM_H

#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

void mpn_printf(mp_limb_t* in, mp_size_t len);
mp_limb_t* mpn_init_cpyz(mpz_t x, size_t n);
mp_limb_t* mpn_init_zero(size_t n);
void mpn2bytes(unsigned char *str, mp_size_t len, const mp_limb_t *s1p, mp_size_t s1n);

#endif