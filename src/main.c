#include <stdio.h>
#include <gmp.h>
#include "ecc.h"
#include "ecc_x.h"
#include "num.h"
#include "dlog.h"

void test_dlog(ecc curve, eccpt G)
{
    eccpt kG;
    ecc_init_pt(kG);

    mpz_t k;
    mpz_t n;
    mpz_init_set_str(k, "2023", 10);
    mpz_init_set_str(n, "857765763956341", 10);

    ecc_mul(curve, kG, G, k);

    unsigned int n_threads = 4;
    dlog(curve, k, G, kG, n, n_threads);
    
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