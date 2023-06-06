#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <gmp.h>
#include "ecc.h"
#include "ecc_x.h"
#include "ex_mpn.h"
#include "dlog.h"
#include "const.h"
#include "mem.h"

const unsigned int BLOCK_BUFFER_SIZE = 1024;
int str_init_readline(char** pbuffer)
{
    char c = 's';
    int buffer_size = BLOCK_BUFFER_SIZE;
    int str_len = 0;
    (*pbuffer) = (char*) malloc_exit_when_null(buffer_size);

    ssize_t nbytes_read;
    while (c != '\0') {
        // Read character
        nbytes_read = read(STDIN_FILENO, &c, 1);
        if (!nbytes_read || c == '\n' || c == '\0')
            c = '\0';
        else if (c < '0' || c > '9')
            continue;

        // Reallocate if needed
        if (str_len == buffer_size) {
            buffer_size += BLOCK_BUFFER_SIZE;
            (*pbuffer) = (char*) realloc(*pbuffer, buffer_size);
        }

        (*pbuffer)[str_len++] = c;
    }

    return str_len;
}

unsigned int get_number_of_threads(int argc, char** argv)
{   
    long int NUM_THREADS = DEFAULT_NUM_THREADS;
    if (argc >= 2) {
        NUM_THREADS = strtol(argv[1], NULL, 10);
        if (NUM_THREADS == 0 || NUM_THREADS == LONG_MAX || NUM_THREADS == LONG_MIN)
            NUM_THREADS = DEFAULT_NUM_THREADS;
        else
            NUM_THREADS = abs(NUM_THREADS);
    }
    return (unsigned int) NUM_THREADS;
}

unsigned int get_cache_size(int argc, char** argv)
{   
    // todo: change this.
    long int CACHE_SIZE = DEFAULT_NUM_THREADS;
    if (argc >= 2) {
        CACHE_SIZE = strtol(argv[2], NULL, 10);
        if (CACHE_SIZE == 0 || CACHE_SIZE == LONG_MAX || CACHE_SIZE == LONG_MIN)
            CACHE_SIZE = DEFAULT_CACHE_SIZE;
        else
            CACHE_SIZE = abs(CACHE_SIZE);
    }
    return (unsigned int) CACHE_SIZE;
}

unsigned int get_random_size(int argc, char** argv)
{   
    long int RAND_SIZE = DEFAULT_NRANDPOINTS;
    if (argc >= 2) {
        RAND_SIZE = strtol(argv[3], NULL, 10);
        if (RAND_SIZE == 0 || RAND_SIZE == LONG_MAX || RAND_SIZE == LONG_MIN)
            RAND_SIZE = DEFAULT_NUM_THREADS;
        else
            RAND_SIZE = abs(RAND_SIZE);
    }
    return (unsigned int) RAND_SIZE;
}

void main(int argc, char** argv)
{
    // Get value from arguments.
    unsigned int NUM_THREADS = get_number_of_threads(argc, argv);
    unsigned int CACHE_SIZE  = get_cache_size(argc, argv);
    unsigned int RAND_SIZE   = get_random_size(argc, argv);

    // Parse curve values from stdin.
    char* curve_a;
    char* curve_b;
    char* curve_p;
    char* Gx;
    char* Gy;
    char* kGx;
    char* kGy;
    char* n_str;

    str_init_readline(&curve_a);
    str_init_readline(&curve_b);
    str_init_readline(&curve_p);
    str_init_readline(&Gx);
    str_init_readline(&Gy);
    str_init_readline(&kGx);
    str_init_readline(&kGy);
    str_init_readline(&n_str);

    // parse curve.
    ecc curve;
    ecc_init(
        curve, 
        curve_a, // a 
        curve_b, // b
        curve_p  // p
    );

    // G
    eccpt G;
    ecc_init_pt_str(
        curve, G,
        Gx,   // x
        Gy,   // y
        NULL  // z
    );
    
    // k*G
    eccpt kG;
    ecc_init_pt_str(
        curve, kG,
        kGx,
        kGy,
        NULL
    );

    // order of G.
    mpz_t n;
    mpz_init_set_str(n, n_str, 10);

    // dlog() start.
    mpz_t k;
    mpz_init(k);
    if (dlog2(
            curve, 
            
            k, 
            G, kG, 
            n, 
            
            NUM_THREADS, 
            CACHE_SIZE, 
            RAND_SIZE
        ) == DLOG_SUCCESS
    ) 
    {
        mpz_out_str(stdout, 10, k);
        printf("\n");
    }
    else
        printf("None\n");

    // cleanup.
    mpz_clear(k);
    mpz_clear(n);

    ecc_free_pt(G);
    ecc_free_pt(kG);
    ecc_free(curve);

    free(curve_a);
    free(curve_b);
    free(curve_p);
    free(Gx);
    free(Gy);
    free(kGx);
    free(kGy);
    free(n_str);
}