#include <stdio.h>
#include <gmp.h>
#include "ecc.h"
#include "ecc_proj.h"
#include "ex_mpn.h"
#include "ex_mpz.h"
#include "dlog.h"

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

void test1()
{
    ecc curve;
    ecc_init(
        curve, 
        "435886898129897736",    // a 
        "255399037257419668",    // b
        "493424529467329699"     // p
    );

    printf("***************** SAGEMATH DEBUG CODE *****************\n");
    printf("a = "); mpz_out_str(stdout, 10, curve->a); printf("\n");
    printf("b = "); mpz_out_str(stdout, 10, curve->b); printf("\n");
    printf("p = "); mpz_out_str(stdout, 10, curve->p); printf("\n");
    printf("E = EllipticCurve(GF(p), [a, b])\n");

    mpz_t n;
    mpz_init_set_str(n, "493424528226388063", 10);
    printf("n = ");
    mpz_out_str(stdout, 10, n);
    printf("\n");

    eccpt G;
    ecc_init_pt(G);
    ecc_random_pt(curve, G);
    printf("G = ");
    eccpt_sage_printf(G);

    mpz_t k;
    mpz_init(k);
    mpz_dev_urandomm(k, n);

    eccpt kG;
    ecc_init_pt(kG);
    ecc_mul_noverify(curve, kG, G, k);
    printf("kG = ");
    eccpt_sage_printf(kG);

    printf("EE = pari.ellinit([E.a4(), E.a6()], p)\n");
    printf("pari.elllog(EE, kG, G, n)\n");

    printf("*******************************************************\n\n");

    dlog(
        curve, 
        k,
        G, kG, 
        n,
        4,
        5,
        20
    );

    printf("k = ");
    mpz_out_str(stdout, 10, k);
    printf("\n");

    mpz_clear(k);
    mpz_clear(n);
    ecc_free_pt(G);
    ecc_free_pt(kG);
    ecc_free(curve);
}

int main()
{
    // test4();
    // test5();
    // test6();
    // test7();
    // test8();
    // test9();
    // test10();
    test1();
}