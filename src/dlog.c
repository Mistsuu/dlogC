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
    assertf(ecc_verify_pt(curve, G), "[" SHARED_LIB_NAME "] ERR: Please provide point G that is on the curve.");
    assertf(ecc_verify_pt(curve, kG), "[" SHARED_LIB_NAME "] ERR: Please provide point kG that is on the curve.");

    // Is order positive and not zero?
    assertf(mpz_cmp_si(G_mult_order, 0) > 0, "[" SHARED_LIB_NAME "] ERR: Please provide a positive value G_mult_order to compute dlog!");

    // Checks G * order = O
    // Checks kG * order = O
    eccpt O;
    ecc_init_pt(O);
    ecc_mul_noverify(curve, O, G, G_mult_order);
    assertf(mpz_cmp_ui(O->z, 0) == 0, "[" SHARED_LIB_NAME "] ERR: Please provide G_mult_order that G_mult_order * G = infinity to compute dlog!");
    ecc_mul_noverify(curve, O, kG, G_mult_order);
    assertf(mpz_cmp_ui(O->z, 0) == 0, "[" SHARED_LIB_NAME "] ERR: Please provide G_mult_order that G_mult_order * kG = infinity to compute dlog!");

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