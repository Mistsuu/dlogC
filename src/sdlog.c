#include "sdlog.h"
#include "dlog.h"
#include "ex_mpn.h"
#include "const.h"
#include "mem.h"

void user_interrupt_handler(
    int signum
)
{
    printf("[" SHARED_LIB_NAME "] Caught SIGINT (signum = %d)! Exiting in peace...\n", signum);
    exit(-1);
}

int sdlog(
    // Curve parameters
    char* str_curve_a,
    char* str_curve_b,
    char* str_curve_p,

    // To be overwritten
    char** pstr_k,

    char* str_Gx,
    char* str_Gy,
    char* str_Gz,
    char* str_kGx,
    char* str_kGy,
    char* str_kGz,
    char* str_n,

    // Configs
    unsigned long n_threads,
    unsigned long alpha,
    unsigned long n_rand_items
)
{
    // Register user's interrupt :)
    signal(SIGINT, user_interrupt_handler);

    ecc curve;
    ecc_init(
        curve, 
        str_curve_a, // a 
        str_curve_b, // b
        str_curve_p  // p
    );

    mpz_t k;
    mpz_init(k);

    eccpt G;
    ecc_init_pt_str(
        curve, G,
        str_Gx,   // x
        str_Gy,   // y
        str_Gz    // z
    );

    eccpt kG;
    ecc_init_pt_str(
        curve, kG,
        str_kGx,
        str_kGy,
        str_kGz
    );

    mpz_t n;
    mpz_init_set_str(n, str_n, 10);

    // Calling the inner dlog().
    int dlog_success = (dlog(
                            curve,
                            k,
                            G, kG,
                            n,
                            
                            n_threads,
                            alpha == 0 
                                ? (unsigned long)mpz_sizeinbase(n, 2) * 3
                                : alpha,
                            n_rand_items
                        ) == DLOG_SUCCESS); 

    // Returns "None" if cannot find.
    if (!dlog_success) {
        *pstr_k = (char*) malloc_exit_when_null(sizeof(char) * 5);
        (*pstr_k)[0] = 'N';
        (*pstr_k)[1] = 'o';
        (*pstr_k)[2] = 'n';
        (*pstr_k)[3] = 'e';
        (*pstr_k)[4] = '\0';
        goto sdlog_cleanup;
    }

    // Convert k to base-10 string
    size_t str_k_len = mpz_sizeinbase(k, 10) + 2;
    *pstr_k = (char*) malloc_exit_when_null(sizeof(char) * str_k_len);
    mpz_get_str(*pstr_k, 10, k);

    // Cleanup.
sdlog_cleanup:
    mpz_clear(k);
    mpz_clear(n);

    ecc_free_pt(G);
    ecc_free_pt(kG);
    ecc_free(curve);

    return dlog_success;
}

void sdlog_free(
    char* str_k
)
{
    free(str_k);
}