#ifndef RAND_H
#define RAND_H

#include <gmp.h>
#include <stdio.h>

void mpz_dev_urandomm(mpz_t rop, const mpz_t n);

#endif