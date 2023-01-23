#ifndef DLOG_H
#define DLOG_H

#include <gmp.h>
#include <assert.h>
#include <stdint.h>
#include "ecc.h"
#include "ecc_x.h"

char* dlog_alloc_buffer(mpz_t nitems, size_t index_size, size_t item_size);
int dlog(ecc curve, mpz_t k, eccpt G, eccpt kG, mpz_t upper_k, unsigned int n_threads);

#endif