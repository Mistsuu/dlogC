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

#define eccpt_sage_printf(P)                      \
do {                                              \
    printf("E(");                                 \
    mpz_out_str(stdout, 10, P->x); printf(", ");  \
    mpz_out_str(stdout, 10, P->y); printf(", ");  \
    mpz_out_str(stdout, 10, P->z); printf(")\n"); \
} while (0)

#define mpnpt_sage_printf(Px, Py, Pz)   \
do {                                    \
    printf("E(");                       \
    mpn_printf(Px, n); printf(", ");    \
    mpn_printf(Py, n); printf(", ");    \
    mpn_printf(Pz, n); printf(")\n");   \
} while (0)

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
    // Code to generate curve in Sagemath

    printf("***************** SAGEMATH DEBUG CODE *****************\n");
    printf("a = "); mpz_out_str(stdout, 10, curve->a); printf("\n");
    printf("b = "); mpz_out_str(stdout, 10, curve->b); printf("\n");
    printf("p = "); mpz_out_str(stdout, 10, curve->p); printf("\n");
    printf("E = EllipticCurve(GF(p), [a, b])\n");
    printf("*******************************************************\n\n");

    n = mpz_size(curve->p);

    mpz_t R;
    mpz_t P;   // = -p^-1 mod R
    mpz_t b3R; // = 3*b
    mpz_t aR;
    mpz_init(R);
    mpz_init(P);
    mpz_init(aR);
    mpz_init_set(b3R, curve->b);
    mpz_init_set(aR, curve->a);
    mpz_set_ui(R, 1);
    mpz_mul_2exp(R, R, n*mp_bits_per_limb);
    mpz_invert(P, curve->p, R);
    mpz_sub(P, R, P);
    mpz_mul(aR, aR, R);
    mpz_mod(aR, aR, curve->p);
    mpz_mul_ui(b3R, b3R, 3);
    mpz_mul(b3R, b3R, R);
    mpz_mod(b3R, b3R, curve->p);

    eccpt G;
    ecc_init_pt(G);
    ecc_random_pt(curve, G);
    printf("G = ");
    eccpt_sage_printf(G);

    eccpt S;
    ecc_init_pt(S);
    ecc_random_pt(curve, S);
    printf("S = ");
    eccpt_sage_printf(S);

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
    eccpt_sage_printf(kG);
    
    eccpt GG;
    ecc_init_pt(GG);
    ecc_add_noverify(curve, GG, G, G);
    printf("2G = ");
    eccpt_sage_printf(GG);

    eccpt STmtp;
    ecc_init_pt(STmtp); 
    ecc_add_noverify(curve, STmtp, G, S);
    printf("G + S = ");
    eccpt_sage_printf(STmtp);

    mp_limb_t* Rx = mpn_init_zero(n);
    mp_limb_t* Ry = mpn_init_zero(n);
    mp_limb_t* Rz = mpn_init_zero(n);
    mp_limb_t* Px = mpn_init_cpyz(G->x, n);
    mp_limb_t* Py = mpn_init_cpyz(G->y, n);
    mp_limb_t* Pz = mpn_init_cpyz(G->z, n);
    mp_limb_t* Qx = mpn_init_cpyz(S->x, n);
    mp_limb_t* Qy = mpn_init_cpyz(S->y, n);
    mp_limb_t* Qz = mpn_init_cpyz(S->z, n);
    mp_limb_t* k  = mpn_init_cpyz(k_, n);
    mp_limb_t* curve_aR  = mpn_init_cpyz(aR, n);
    mp_limb_t* curve_b3R = mpn_init_cpyz(b3R, n);
    mp_limb_t* curve_p   = mpn_init_cpyz(curve->p, n);
    mp_limb_t* curve_P   = mpn_init_cpyz(P, n);

    ecc_ptemp T;
    ecc_init_ptemp(T, n);

    ecc_pmul(
        Rx, Ry, Rz,
        Px, Py, Pz,
        k,
        curve_aR,
        curve_b3R,
        curve_p,
        curve_P,
        n,
        T
    );

    // ecc_padd(
    //     Rx, Ry, Rz,
    //     Px, Py, Pz,
    //     Qx, Qy, Qz,
    //     curve_aR,
    //     curve_b3R,
    //     curve_p,
    //     curve_P,
    //     n,
    //     T
    // );

    // ecc_pdbl(
    //     Rx, Ry, Rz,
    //     Px, Py, Pz,
    //     curve_aR,
    //     curve_b3R,
    //     curve_p,
    //     curve_P,
    //     n,
    //     T
    // );

    mpnpt_sage_printf(Rx, Ry, Rz);
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