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
    // ----------------------------- G, kG on curve? -----------------------------
    if (!ecc_verify_pt(curve, G)) {
        #ifdef DLOG_VERBOSE
            printf("[debug] G is not on the curve. Exiting...\n");
        #endif
        return DLOG_POINT_NOT_ON_CURVE;
    }

    if (!ecc_verify_pt(curve, kG)) {
        #ifdef DLOG_VERBOSE
            printf("[debug] kG is not on the curve. Exiting...\n");
        #endif
        return DLOG_POINT_NOT_ON_CURVE;
    }

    // ---------------------- is order positive and prime? -----------------------
    if (mpz_cmp_si(G_mult_order, 0) <= 0) {
        #ifdef DLOG_VERBOSE
            printf("[debug] G_mult_order is negative. Exiting...\n");
        #endif
        return DLOG_FAULTY_POINT_ORDER;
    }

    if (!mpz_probab_prime_p(G_mult_order, BASESIZE_PRIME_CHECKER)) {
        #ifdef DLOG_VERBOSE
            printf("[debug] G_mult_order isn't prime. Exiting...\n");
        #endif
        return DLOG_FAULTY_POINT_ORDER;
    }

    // -------------------------- checks G * order = O ---------------------------
    eccpt O;
    ecc_init_pt(O);
    ecc_mul_noverify(curve, O, G, G_mult_order);
    if (mpz_cmp_ui(O->z, 0) != 0) {
        #ifdef DLOG_VERBOSE
            printf("[debug] G_mult_order * G != infinity. Exiting...\n");
        #endif
        ecc_free_pt(O);
        return DLOG_FAULTY_POINT_ORDER;
    }

    ecc_free_pt(O);
    return DLOG_MOVE_TO_NEXT_STEP;
}

int dlog_fast_solve_if_possible(
    ecc curve,
    mpz_t k, 
    eccpt G, eccpt kG,
    mpz_t G_mult_order
)
{
    // -------------------------------------------------------------------------------------
    //      Case 0: kG = O
    // -------------------------------------------------------------------------------------
    if (mpz_cmp_ui(kG->z, 0) == 0) {
            #ifdef DLOG_VERBOSE
            printf("[debug] kG == O. Automatically set k = 0.\n");
        #endif
        mpz_set_ui(k, 0);
        return DLOG_SUCCESS;
    }

    // -------------------------------------------------------------------------------------
    //      Case 1: G == O and kG != O
    // -------------------------------------------------------------------------------------
    if (mpz_cmp_ui(G->z, 0) == 0) {
        #ifdef DLOG_VERBOSE
            printf("[debug] G == O and kG != O. No solution exists.\n");
        #endif
        return DLOG_NOT_FOUND_DLOG;
    }

    // -------------------------------------------------------------------------------------
    //      Case 2: kG.order() > G.order()
    // -------------------------------------------------------------------------------------
    eccpt O;
    ecc_init_pt(O);
    ecc_mul_noverify(curve, O, kG, G_mult_order);
    if (mpz_cmp_ui(O->z, 0) != 0) {
        #ifdef DLOG_VERBOSE
            printf("[debug] n*G == infinity but n*kG != infinity. No solution exists.\n");
        #endif
        ecc_free_pt(O);
        return DLOG_NOT_FOUND_DLOG;
    }

    // -------------------------------------------------------------------------------------
    //      Case 3: Using Weil Pairing to 
    //              check if there's solution?
    // -------------------------------------------------------------------------------------
    mpz_t E;
    mpz_init(E);
    ecc_weil_pairing_noverify(curve, E, G, kG, G_mult_order);
    if (mpz_cmp_ui(E, 1) != 0) {
        #ifdef DLOG_VERBOSE
            printf("[debug] weil_pairing(G, kG) != 1. kG is not a multiple of G.\n");
        #endif
        mpz_clear(E);
        ecc_free_pt(O);
        return DLOG_NOT_FOUND_DLOG;
    }

    mpz_clear(E);
    ecc_free_pt(O);
    return DLOG_MOVE_TO_NEXT_STEP;
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
    #ifdef DLOG_VERBOSE
        printf("[debug] curve: \n");
        printf("[debug]    ");
        ecc_printf(curve);
        printf("\n");
        printf("[debug] G: \n");
        printf("[debug]    ");
        ecc_printf_pt(G);
        printf("\n");
        printf("[debug] kG: \n");
        printf("[debug]    ");
        ecc_printf_pt(kG);
        printf("\n");
        printf("[debug] G_mult_order = ");
        mpz_out_str(stdout, 10, G_mult_order);
        printf("\n");
        printf("[debug] n_threads = %d\n", n_threads);
        printf("[debug] cache_size = %d\n", cache_size);
    #endif

    int dlog_status;

    // -------------------------------------------------------------------------------------
    //      Verify if input is good.
    // -------------------------------------------------------------------------------------
    dlog_status = dlog_validate_input(curve, G, kG, G_mult_order);
    if (dlog_status != DLOG_MOVE_TO_NEXT_STEP)
        return dlog_status;

    // -------------------------------------------------------------------------------------
    //      Quickly determine if we can/or can't solve
    //      and provide an answer if possible. The reason
    //      for this is that the Pollard-rho algorithm
    //      is probabilistic rather than deterministic. Doing
    //      this will prevent the algorithm from spiraling
    //      into an infinite loop.
    // -------------------------------------------------------------------------------------
    dlog_status = dlog_fast_solve_if_possible(curve, k, G, kG, G_mult_order);
    if (dlog_status != DLOG_MOVE_TO_NEXT_STEP)
        return dlog_status;
}