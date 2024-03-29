#include "ex_mpz.h"
#include "mem.h"

// =================================================================================
//                                 MEMORY STUFFS
// =================================================================================

#define mpz_size_bytes(n) mpz_sizeinbase(n, 256)
void mpz_dev_urandomm(mpz_t rop, const mpz_t n)
{
    FILE *fptr = fopen("/dev/urandom", "r");
    if (!fptr) {
        printf("[error] Error! Cannot read from random source!\n");
        exit(-1);
    }

    // Read random bytes
    // to mp_limb_t.
    size_t item_size_limbs = mpz_size(n);
    mp_limb_t* tmp_limbs; 
    tmp_limbs = (mp_limb_t*) malloc_exit_when_null(item_size_limbs * sizeof(mp_limb_t));
    if (fread(tmp_limbs, sizeof(mp_limb_t), item_size_limbs, fptr) != item_size_limbs) {
        printf("[error] Error! Memory source is insufficient!\n");
        exit(-1);
    }

    // Convert to mpz_t
    mpz_t tmp;
    mpz_set(rop, mpz_roinit_n(tmp, tmp_limbs, item_size_limbs));
    mpz_mod(rop, rop, n);

    // Free memory.
    free(tmp_limbs);
    fclose(fptr);
}

void mpz_init_set_mpn(mpz_t rop, mp_limb_t* xp, mp_size_t n)
{
    mpz_init(rop);
    mp_limb_t* tp = mpz_limbs_modify(rop, n);
    mpn_copyd(tp, xp, n);
    while (n && !xp[n-1])
        n--;
    mpz_limbs_finish(rop, n);
}

void mpz_set_mpn(mpz_t rop, mp_limb_t* xp, mp_size_t n)
{
    mp_limb_t* tp = mpz_limbs_modify(rop, n);
    mpn_copyd(tp, xp, n);
    while (n && !xp[n-1])
        n--;
    mpz_limbs_finish(rop, n);
}

// =================================================================================
//                                 ARITHMETICS STUFFS
// =================================================================================

/* Source code taken from https://gmplib.org/list-archives/gmp-discuss/2007-May/002772.html */

/* mpz_sqrtm -- modular square roots using Shanks-Tonelli

Copyright 2006 Free Software Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
for more details.

You should have received a copy of the GNU Lesser General Public License along
with the GNU MP Library; see the file COPYING.LIB.  If not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
USA. */

/* Solve the modular equatioon x^2 = n (mod p) using the Shanks-Tonelli
 * algorihm. x will be placed in q and 1 returned if the algorithm is
 * successful. Otherwise 0 is returned (currently in case n is not a quadratic
 * residue mod p). A check is done if p = 3 (mod 4), in which case the root is
 * calculated as n ^ ((p+1) / 4) (mod p).
 *
 * Note that currently mpz_legendre is called to make sure that n really is a
 * quadratic residue. The check can be skipped, at the price of going into an
 * eternal loop if called with a non-residue.
 */
int mpz_sqrtm(mpz_t q, const mpz_t n, const mpz_t p) {
    mpz_t w, n_inv, y;
    unsigned int i, s;

    if(mpz_divisible_p(n, p)) {         /* Is n a multiple of p?            */
        mpz_set_ui(q, 0);               /* Yes, then the square root is 0.  */
        return 1;                       /* (special case, not caught        */
    }                                   /* otherwise)                       */
    if(mpz_legendre(n, p) != 1)         /* Not a quadratic residue?         */
        return 0;                       /* No, so return error              */
    if(mpz_tstbit(p, 1) == 1) {         /* p = 3 (mod 4) ?                  */
        mpz_set(q, p);
        mpz_add_ui(q, q, 1);
        mpz_fdiv_q_2exp(q, q, 2);
        mpz_powm(q, n, q, p);           /* q = n ^ ((p+1) / 4) (mod p)      */
        return 1;
    }
    mpz_init(y);
    mpz_init(w);
    mpz_init(n_inv);
    mpz_set(q, p);
    mpz_sub_ui(q, q, 1);                /* q = p-1                          */
    s = 0;                              /* Factor out 2^s from q            */
    while(mpz_tstbit(q, s) == 0) s++;
    mpz_fdiv_q_2exp(q, q, s);           /* q = q / 2^s                      */
    mpz_set_ui(w, 2);                   /* Search for a non-residue mod p   */
    while(mpz_legendre(w, p) != -1)     /* by picking the first w such that */
        mpz_add_ui(w, w, 1);            /* (w/p) is -1                      */
    mpz_powm(w, w, q, p);               /* w = w^q (mod p)                  */
    mpz_add_ui(q, q, 1);
    mpz_fdiv_q_2exp(q, q, 1);           /* q = (q+1) / 2                    */
    mpz_powm(q, n, q, p);               /* q = n^q (mod p)                  */
    mpz_invert(n_inv, n, p);
    for(;;) {
        mpz_powm_ui(y, q, 2, p);        /* y = q^2 (mod p)                  */
        mpz_mul(y, y, n_inv);
        mpz_mod(y, y, p);               /* y = y * n^-1 (mod p)             */
        i = 0;
        while(mpz_cmp_ui(y, 1) != 0) {
            i++;
            mpz_powm_ui(y, y, 2, p);    /*  y = y ^ 2 (mod p)               */
        }
        if(i == 0) {                    /* q^2 * n^-1 = 1 (mod p), return   */
            mpz_clear(y);
            mpz_clear(w);
            mpz_clear(n_inv);
            return 1;
        }
        if(s-i == 1) {                  /* In case the exponent to w is 1,  */
            mpz_mul(q, q, w);           /* Don't bother exponentiating      */
        } else {
            mpz_powm_ui(y, w, 1 << (s-i-1), p);
            mpz_mul(q, q, y);
        }
        mpz_mod(q, q, p);               /* r = r * w^(2^(s-i-1)) (mod p)    */
    }
}