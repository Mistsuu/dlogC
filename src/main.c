#include <stdio.h>
#include <gmp.h>
#include "ecc.h"
#include "ecc_x.h"
#include "num.h"
#include "dlog.h"

void test3()
{
    mpz_t p;
    size_t n;
    mpz_init(p);
    mpz_set_str(p, "21272169233168579221109700210286068195545041594730488403677158434097465979633", 10);
    n = mpz_size(p);

    mp_limb_t* _p = mpn_init_cpyz(p, n);
    mp_limb_t* _x = mpn_init_zero(n);
    mp_limb_t* _y = mpn_init_zero(n);
    mp_limb_t* _z = mpn_init_zero(2*n);
    mp_limb_t* _T = mpn_init_zero(2*n);
    mp_limb_t* _A = mpn_init_zero(2*n);
    mp_limb_t* _B = mpn_init_zero(2*n);
    mp_size_t t;

    mpn_random(_x, n);
    mpn_random(_y, n);
    mpn_mul_n(_z, _x, _y, n);
    
    for (int _ = 0; _ < 10000000; ++_)
    {
        mpn_tdiv_qr(
            _z, _T, 0,
            _z, 2*n,
            _p, n
        );

        mpn_gcdext( 
            _A, _B, &t,
            _x, n,
            _p, n
        );
    }

    free(_p);
    mpz_clear(p);
}

void test4()
{
    ecc curve;
    eccpt G;
    ecc_init(
        curve, 
        "1986076773346522069111732327339",    // a 
        "808177731529494834911895879646",     // b
        "13276420418771432419898581447951"    // p
    );

    ecc_init_pt(G);
    ecc_random_pt(curve, G);
    ecc_printf_pt(G);
    printf("\n");

    ecc_random_pt(curve, G);
    ecc_printf_pt(G);
    printf("\n");

    ecc_random_pt(curve, G);
    ecc_printf_pt(G);
    printf("\n");

    ecc_free_pt(G);
    ecc_free(curve);
}

int main()
{
    // test3();
    test4();
}