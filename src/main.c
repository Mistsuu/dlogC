#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <gmp.h>
#include "ecc.h"
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

unsigned long parse_arg_ulong(char* arg, unsigned long default_val)
{
    unsigned long val = default_val;
    if (arg) {
        val = strtoul(arg, NULL, 10);
        if (val == 0 || val == ULONG_MAX)
            val = default_val;
    }
    return val;
}

int main(int argc, char** argv)
{
    // ----------------------------- parse curve values from stdin. -----------------------------
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

    // ----------------------------- parse value from argc, argv. -----------------------------
    unsigned long NUM_THREADS     = DEFAULT_NUM_THREADS;
    unsigned long NUM_RAND_ITEMS  = DEFAULT_NUM_RAND_ITEMS;
    unsigned long ALPHA           = (unsigned long)mpz_sizeinbase(n, 2) * 3;
    
    int opt;
    while ((opt = getopt(argc, argv, "t:a:r:h")) != -1) {
        switch (opt) {
            case 't':
                NUM_THREADS = parse_arg_ulong(optarg, DEFAULT_NUM_THREADS);
                break;
            case 'a':
                ALPHA = parse_arg_ulong(optarg, (unsigned long)mpz_sizeinbase(n, 2) * 3);
                break;
            case 'r':
                NUM_RAND_ITEMS = parse_arg_ulong(optarg, DEFAULT_NUM_RAND_ITEMS);
                break;
            default:
                fprintf(stderr, "[usage]: %s [-t num_threads=4] [-a alpha=3*log2(n)] [-r nrandpoints=20]\n", argv[0]);
                exit(-1);
        }
    }
    
    // ---------------------------------- run dlog() ----------------------------------------
    mpz_t k;
    mpz_init(k);
    if (dlog(
            curve, 
            k, 
            G, kG, 
            n, 
            NUM_THREADS, 
            ALPHA, 
            NUM_RAND_ITEMS
        ) == DLOG_SUCCESS
    ) 
    {
        mpz_out_str(stdout, 10, k);
        printf("\n");
    }
    else
        printf("None\n");

    // ---------------------------------- cleaning... ----------------------------------------
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

    return 0;
}