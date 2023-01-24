#ifndef NUM_H
#define NUM_H

#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

void mpn_printf(mp_limb_t* in, mp_size_t len);
mp_limb_t* mpz_limbs_init_cpy(mpz_t x, size_t n);
mp_limb_t* mpz_limbs_init_zero(size_t n);

#endif