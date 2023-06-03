#include <stdio.h>
#include <gmp.h>
#include "ecc.h"
#include "ecc_x.h"
#include "ex_mpn.h"
#include "olddlog.h"
#include "const.h"

void test_dlog(ecc curve, eccpt G)
{
    eccpt kG;
    ecc_init_pt(kG);

    mpz_t k;
    mpz_t n;
    mpz_init_set_str(k, "22535525235", 10);
    mpz_init_set_str(n, "857765763956341", 10);

    ecc_mul(curve, kG, G, k);

    printf("[i] Finding k = dlog(G, kG) for point:\n");
    printf("[i]    G = "); ecc_printf_pt(G); printf("\n");
    printf("[i]    kG = "); ecc_printf_pt(kG); printf("\n");
    printf("[i] on curve: "); printf("\n");
    printf("[i]    "); ecc_printf(curve); printf("\n");

    unsigned int n_threads = 4;
    if (dlog(curve, k, G, kG, n, n_threads) == DLOG_SUCCESS)
    {
        printf("[i] k = "); 
        mpz_out_str(stdout, 10, k);
        printf("\n");
    }
    else
    {
        printf("[i] Cannot find dlog!\n");
    }
    
    mpz_clear(k);
    mpz_clear(n);
    ecc_free_pt(kG);
}

void main()
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

    test_dlog(curve, G);

    ecc_free_pt(G);
    ecc_free(curve);
}