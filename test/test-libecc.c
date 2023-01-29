#include <stdio.h>
#include <gmp.h>
#include "ecc.h"
#include "ecc_x.h"
#include "num.h"
#include "dlog.h"

void test1(ecc curve, eccpt G)
{
    eccpt R;
    mpz_t k;
    mpz_init(k);
    ecc_init_pt(R);

    mpz_set_str(k, "2023", 10);
    ecc_mul(curve, R, G, k);

    ecc_printf_pt(G);
    printf("\n");
    ecc_printf_pt(R);
    printf("\n");

    mpz_clear(k);
    ecc_free_pt(R);
}

void test2(ecc curve, eccpt G)
{
    mp_size_t n = curve->p->_mp_size;
    mp_limb_t k[2] = {2023, 0};

    mp_limb_t* Gx = (mp_limb_t*) malloc(sizeof(mp_limb_t) * n);
    mp_limb_t* Gz = (mp_limb_t*) malloc(sizeof(mp_limb_t) * n);
    mp_limb_t* Rx = (mp_limb_t*) malloc(sizeof(mp_limb_t) * n);
    mp_limb_t* Rz = (mp_limb_t*) malloc(sizeof(mp_limb_t) * n);
    mp_limb_t* curve_p = (mp_limb_t*) malloc(sizeof(mp_limb_t) * n);
    mp_limb_t* curve_a = (mp_limb_t*) malloc(sizeof(mp_limb_t) * n);
    mp_limb_t* curve_b = (mp_limb_t*) malloc(sizeof(mp_limb_t) * n);

    mpn_zero(Gx, n);
    mpn_zero(Gz, n);
    mpn_zero(curve_p, n);
    mpn_zero(curve_a, n);
    mpn_zero(curve_b, n);
    mpn_copyd(Gx, G->x->_mp_d, G->x->_mp_size);
    mpn_copyd(Gz, G->z->_mp_d, G->z->_mp_size);
    mpn_copyd(curve_p, curve->p->_mp_d, curve->p->_mp_size);
    mpn_copyd(curve_a, curve->a->_mp_d, curve->a->_mp_size);
    mpn_copyd(curve_b, curve->b->_mp_d, curve->b->_mp_size);

    ecc_xtemp T;
    ecc_init_xtemp(T, n);

    for (int _ = 0; _ < 1000000; ++_)
        ecc_xmul(
            Rx, Rz,
            Gx, Gz,
            k,
            curve_a,
            curve_b,
            curve_p,
            n,
            T
        );

    printf("GF(p)("); mpn_printf(Rx, n); printf(")/"); mpn_printf(Rz, n); printf("\n");

    ecc_free_xtemp(T);
    free(Rx);
    free(Rz);
    free(Gx);
    free(Gz);
    free(curve_p);
    free(curve_a);
    free(curve_b);
}

int main()
{
    ecc curve;
    eccpt G;
    ecc_init(
        curve, 
        "1986076773346522069111732327339",    // a 
        "808177731529494834911895879646",     // b
        "13276420418771432419898581447951"    // p
    );

    ecc_init_pt_str(
        curve, G,
        "12752653901711390718579996242468",   // x
        "9102988295173351464328400869432",    // y
        NULL                                  // z
    );

    test1(curve, G);
    test2(curve, G);

    ecc_free_pt(G);
    ecc_free(curve);

}