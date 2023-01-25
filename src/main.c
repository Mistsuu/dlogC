#include <stdio.h>
#include <gmp.h>
#include "ecc.h"
#include "ecc_x.h"
#include "num.h"
#include "dlog.h"
#include "const.h"

const unsigned int NUM_THREADS = 4;

void main()
{
    ecc curve;
    ecc_init(
        curve, 
        "1986076773346522069111732327339",    // a 
        "808177731529494834911895879646",     // b
        "13276420418771432419898581447951"    // p
    );

    eccpt G;
    ecc_init_pt_str(
        curve, G,
        "12752653901711390718579996242468",   // x
        "9102988295173351464328400869432",    // y
        NULL                                  // z
    );

    eccpt kG;
    ecc_init_pt_str(
        curve, kG,
        "160854798263565084664403423288",
        "2332898824679189780448318708917",
        NULL
    );

    printf("[i] Finding k = dlog(G, kG) for point:\n");
    printf("[i]    G = "); ecc_printf_pt(G); printf("\n");
    printf("[i]    kG = "); ecc_printf_pt(kG); printf("\n");
    printf("[i] on curve: "); printf("\n");
    printf("[i]    "); ecc_printf(curve); printf("\n");

    // order of G.
    mpz_t n;
    mpz_init_set_str(n, "857765763956341", 10);

    mpz_t k;
    mpz_init(k);
    if (dlog(curve, k, G, kG, n, NUM_THREADS) == DLOG_SUCCESS)
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
    ecc_free_pt(G);
    ecc_free_pt(kG);
    ecc_free(curve);
}