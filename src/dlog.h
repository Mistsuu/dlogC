#ifndef DLOG_H
#define DLOG_H

#include <gmp.h>
#include <assert.h>
#include <stdint.h>
#include "ecc.h"
#include "ecc_x.h"

mp_limb_t* dlog_allocate(mpz_t nitems_alloc, size_t n_size_limbs, size_t p_size_limbs);
int dlog(ecc curve, mpz_t k, eccpt G, eccpt kG, mpz_t n);

#endif