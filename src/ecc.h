#ifndef ECC_H
#define ECC_H

#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <ex_assert.h>

/*
    ecc.h:
        Define elliptic curve and its utilities in the form:
            y^2 = x^3 + ax + b.
        in finite field GF(p).
        
        Using gmp library for multi-precision numbers.
*/

typedef struct ecc_struct
{
    mpz_t a;
    mpz_t b;
    mpz_t p;
} __ecc_struct;

typedef struct eccpt_struct
{
    mpz_t x;
    mpz_t y;
    mpz_t z;
} __eccpt_struct;

typedef __ecc_struct   ecc[1];
typedef __eccpt_struct eccpt[1];

void ecc_init  (ecc curve, const char* a_str, const char* b_str, const char* p_str);
void ecc_printf(ecc curve);
void ecc_free  (ecc curve);

void ecc_init_pt              (eccpt point);
void ecc_init_pt_str          (ecc curve, eccpt point, const char* x_str, const char* y_str, const char* z_str);
void ecc_init_pt_pt           (eccpt dst_point, eccpt src_point);
void ecc_set_pt               (eccpt dst_point, eccpt src_point);
void ecc_set_pt_inf           (eccpt inf_point);
void ecc_normalize_z_pt       (ecc curve, eccpt point);
int  ecc_verify_pt            (ecc curve, eccpt point);
void ecc_random_pt            (ecc curve, eccpt point);
void ecc_printf_pt            (eccpt point);
void ecc_free_pt              (eccpt point);
void ecc_add_noverify         (ecc curve, eccpt pointR, eccpt pointP, eccpt pointQ);
void ecc_mul_noverify         (ecc curve, eccpt pointR, eccpt pointP, mpz_t k);
void ecc_neg_noverify         (ecc curve, eccpt pointR, eccpt pointP);
void ecc_sub_noverify         (ecc curve, eccpt pointR, eccpt pointP, eccpt pointQ);
int  ecc_eq_noverify          (ecc curve, eccpt pointP, eccpt pointQ);
void ecc_add                  (ecc curve, eccpt pointR, eccpt pointP, eccpt pointQ);
void ecc_neg                  (ecc curve, eccpt pointR, eccpt pointP);
void ecc_sub                  (ecc curve, eccpt pointR, eccpt pointP, eccpt pointQ);
void ecc_mul                  (ecc curve, eccpt pointR, eccpt pointP, mpz_t k);
int  ecc_eq                   (ecc curve, eccpt pointP, eccpt pointQ);
void ecc_weil_fP              (ecc curve, mpz_t f, eccpt pointP, eccpt pointR, mpz_t n);
void ecc_weil_gPQ             (ecc curve, mpz_t g, eccpt pointP, eccpt pointQ, eccpt pointR);
void ecc_weil_pairing         (ecc curve, mpz_t E, eccpt pointP, eccpt pointQ, mpz_t n);
void ecc_weil_pairing_noverify(ecc curve, mpz_t E, eccpt pointP, eccpt pointQ, mpz_t n);

#endif