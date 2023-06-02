#include <stdio.h>
#include <gmp.h>
#include "ecc.h"
#include "ecc_x.h"
#include "ex_mpn.h"
#include "dlog.h"
#include "ex_mpz.h"

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
    // mpn_mul_n(_z, _x, _y, n);
    
    for (int _ = 0; _ < 10000000; ++_)
    {
        // mpn_add_n(_z, _x, _y, n);
        mpn_mul_n(_z, _x, _y, n);
        // mpn_tdiv_qr(
        //     _z, _T, 0,
        //     _z, 2*n,
        //     _p, n
        // );

        // mpn_gcdext( 
        //     _A, _B, &t,
        //     _x, n,
        //     _p, n
        // );
    }

    free(_p); 
    free(_x); 
    free(_y); 
    free(_z); 
    free(_T); 
    free(_A); 
    free(_B); 
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

void test5()
{
    mpz_t num;
    mpz_t p;
    mpz_init(num);
    mpz_init_set_str(p, "13276420418771432419898581447951", 10);
    mpz_dev_urandomm(num, p);

    mpz_out_str(stdout, 10, num);
    printf("\n");
    printf("bitlen = %ld\n", mpz_sizeinbase(num, 2));

    mpz_clear(num);
    mpz_clear(p);
}

void test6()
{
    ecc curve;
    eccpt G;
    ecc_init(
        curve, 
        "1986076773346522069111732327339",    // a 
        "808177731529494834911895879646",     // b
        "13276420418771432419898581447951"    // p
    );

    mpz_t n;
    mpz_init_set_str(n, "13276420418771430444004808657717", 10);

    ecc_init_pt(G);
    ecc_random_pt(curve, G);
    ecc_printf_pt(G);
    printf("\n");

    eccpt Z;
    ecc_init_pt(Z);
    ecc_set_pt_inf(Z);

    mpz_t E;
    mpz_init(E);
    ecc_weil_pairing(curve, E, G, Z, n);
    mpz_out_str(stdout, 10, E);
    printf("\n");

    mpz_clear(n);
    mpz_clear(E);
    ecc_free_pt(Z);
    ecc_free_pt(G);
    ecc_free(curve);
}

void test7()
{
    int n = 4;

    // Check what is the percentage of
    // squaring / multiplying.
    mp_limb_t* _x = mpn_init_zero(n);
    mp_limb_t* _y = mpn_init_zero(n);
    mp_limb_t* _z = mpn_init_zero(2*n);
    mpn_random(_x, n);
    mpn_random(_y, n);
    
    for (int _ = 0; _ < 100000000; ++_) {
        // mpn_mul_n(_z, _x, _y, n);
        // mpn_sqr(_z, _x, n);
    }

    free(_x); 
    free(_y); 
    free(_z); 
}

int main()
{
    test3();
    // test4();
    // test5();
    // test6();
    // test7();
}