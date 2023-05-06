#include "sdlog.h"
#include "dlog.h"
#include "num.h"
#include "const.h"
#include "que2.h"
#include "mem.h"

void user_interrupt_handler(
    int signum
)
{
    printf("[libbsgsfp] Caught SIGINT! Exiting in peace...\n");
    exit(-1);
}

int sdlog(
    // Field parameter
    char* str_p,

    // To be modified
    char** pstr_k,

    char* str_G,
    char* str_kG,
    char* str_upper_k,
    
    // Configs
    unsigned int n_threads,
    size_t mem_limit
)
{
    // Register user's interrupt :)
    signal(SIGINT, user_interrupt_handler);

    mpz_t p;
    mpz_init_set_str(p, str_p, 10);
    mpz_t k;
    mpz_init(k);
    mpz_t G;
    mpz_init_set_str(G, str_G, 10);
    mpz_t kG;
    mpz_init_set_str(kG, str_kG, 10);
    mpz_t upper_k;
    mpz_init_set_str(upper_k, str_upper_k, 10);

    // Calling the inner dlog().
    int dlog_success = (dlog(
                            p,
                            k,
                            G,
                            kG,
                            upper_k,
                            n_threads,
                            mem_limit
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
    mpz_clear(p);
    mpz_clear(k);
    mpz_clear(G);
    mpz_clear(kG);
    mpz_clear(upper_k);

    return dlog_success;
}

void sdlog_free(
    char* str_k
)
{
    free(str_k);
}