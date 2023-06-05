#ifndef EXMPZ_H
#define EXMPZ_H

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

void mpz_dev_urandomm(mpz_t rop, const mpz_t n);
void mpz_init_set_mpn(mpz_t rop, mp_limb_t* xp, mp_size_t n);
void mpz_set_mpn     (mpz_t rop, mp_limb_t* xp, mp_size_t n);
int  mpz_sqrtm       (mpz_t q, const mpz_t n, const mpz_t p);

#endif