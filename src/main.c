#include <stdio.h>
#include <gmp.h>
#include "ecc.h"
#include "ecc_x.h"
#include "ecc_proj.h"
#include "ex_mpn.h"
#include "ex_mpz.h"
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

void test8()
{
    mpz_t x;
    mpz_t y;
    mpz_t p;
    mpz_t P; // = -p^-1 mod R
    mpz_t R;
    size_t n;
    mpz_init(x);
    mpz_init(y);
    mpz_init(p);
    mpz_init(P);
    mpz_init(R);

    mpz_set_str(p, "21272169233168579221109700210286068195545041594730488403677158434097465979633", 10);
    n = mpz_size(p);

    mpz_dev_urandomm(x, p);
    mpz_dev_urandomm(y, p);
    mpz_set_ui(R, 1);
    mpz_mul_2exp(R, R, n*mp_bits_per_limb);
    mpz_invert(P, p, R);
    mpz_sub(P, R, P);

    mp_limb_t* xp = mpn_init_cpyz(x, n);
    mp_limb_t* yp = mpn_init_cpyz(y, n);
    mp_limb_t* pp = mpn_init_cpyz(p, n);
    mp_limb_t* Pp = mpn_init_cpyz(P, n);
    mp_limb_t* rp = mpn_init_zero(n);
    mp_limb_t* tp = mpn_init_zero(6*n);
    mp_limb_t* zp = mpn_init_zero(2*n);

    mpn_montgomery_addmod_n(rp, xp, yp, pp, n);

    printf("x = "); mpn_printf(xp, n); printf("\n");
    printf("y = "); mpn_printf(yp, n); printf("\n");
    printf("p = "); mpn_printf(pp, n); printf("\n");
    printf("P = "); mpn_printf(Pp, n); printf("\n");
    printf("r = "); mpn_printf(rp, n); printf("\n");
    printf("R = "); mpz_out_str(stdout, 10, R); printf("\n");
}

void test9()
{
    ecc curve;
    size_t n;
    ecc_init(
        curve, 
        "1986076773346522069111732327339",    // a 
        "808177731529494834911895879646",     // b
        "13276420418771432419898581447951"    // p
    );
    n = mpz_size(curve->p);

    mpz_t R;
    mpz_t P;  // = -p^-1 mod R
    mpz_t b3; // = 3*b
    mpz_init(R);
    mpz_init(P);
    mpz_init_set(b3, curve->b);
    mpz_set_ui(R, 1);
    mpz_mul_2exp(R, R, n*mp_bits_per_limb);
    mpz_invert(P, curve->p, R);
    mpz_sub(P, R, P);
    mpz_mul_ui(b3, b3, 3);
    mpz_mod(b3, b3, curve->p);

    eccpt G;
    ecc_init_pt(G);
    ecc_random_pt(curve, G);
    printf("G = ");
    ecc_printf_pt(G);
    printf("\n");

    mpz_t k_;
    mpz_init(k_);
    mpz_dev_urandomm(k_, curve->p);
    printf("k = ");
    mpz_out_str(stdout, 10, k_);
    printf("\n");

    eccpt kG;
    ecc_init_pt(kG);
    ecc_mul_noverify(curve, kG, G, k_);
    printf("kG = ");
    ecc_printf_pt(kG);
    printf("\n");

    mp_limb_t* Rx = mpn_init_zero(n);
    mp_limb_t* Ry = mpn_init_zero(n);
    mp_limb_t* Rz = mpn_init_zero(n);
    mp_limb_t* Px = mpn_init_cpyz(G->x, n);
    mp_limb_t* Py = mpn_init_cpyz(G->y, n);
    mp_limb_t* Pz = mpn_init_cpyz(G->z, n);
    mp_limb_t* k  = mpn_init_cpyz(k_, n);
    mp_limb_t* curve_a  = mpn_init_cpyz(curve->a, n);
    mp_limb_t* curve_b3 = mpn_init_cpyz(b3, n);
    mp_limb_t* curve_p  = mpn_init_cpyz(curve->p, n);
    mp_limb_t* curve_P  = mpn_init_cpyz(P, n);

    ecc_ptemp T;
    ecc_init_ptemp(T, n);

    ecc_pmul(
        Rx, Ry, Rz,
        Px, Py, Pz,
        k, 
        curve_a,
        curve_b3,
        curve_p,
        curve_P,
        n,
        T
    );

    printf("("); mpn_printf(Rx, n); printf(" : ");
    mpn_printf(Ry, n); printf(" : ");
    mpn_printf(Rz, n); printf(")\n");
}

int main()
{
    // test3();
    // test4();
    // test5();
    // test6();
    // test7();
    // test8();
    test9();
}