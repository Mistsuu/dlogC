#include "dlog.h"
#include "ecc.h"
#include "ecc_proj.h"
#include "ex_mpn.h"
#include "const.h"
#include "que2.h"
#include "mem.h"

int dlog_validate_input(
    ecc curve,
    eccpt G, eccpt kG,
    mpz_t G_mult_order
)
{
    // G, kG on curve?
    assertf(ecc_verify_pt(curve, G), "dlog: G is not on the curve.");
    assertf(ecc_verify_pt(curve, kG), "dlog: kG is not on the curve.");

    // Is order positive and is prime?
    assertf(mpz_cmp_si(G_mult_order, 0) > 0, "dlog: G_mult_order is negative.");
    assertf(mpz_probab_prime_p(G_mult_order, BASESIZE_PRIME_CHECKER), "dlog: G_mult_order isn't prime.");

    // Checks G * order = O
    // Checks kG * order = O
    eccpt O;
    ecc_init_pt(O);
    ecc_mul_noverify(curve, O, G, G_mult_order);
    assertf(mpz_cmp_ui(O->z, 0) == 0, "dlog: G_mult_order * G != infinity.");

    // Free memory
    ecc_free_pt(O);
}

int dlog2(
    ecc curve, 
    mpz_t k, 
    eccpt G, eccpt kG, 
    mpz_t G_mult_order, 

    unsigned int n_threads,
    unsigned int cache_size
)
{

}