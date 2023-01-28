#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <gmp.h>
#include "ecc.h"
#include "ecc_x.h"
#include "num.h"
#include "dlog.h"
#include "const.h"

const unsigned int NUM_THREADS = 4;
const unsigned int BLOCK_BUFFER_SIZE = 10;
int str_init_readline(char** pbuffer)
{
    char c = 's';
    int buffer_size = BLOCK_BUFFER_SIZE;
    int str_len = 0;
    (*pbuffer) = (char*) malloc(buffer_size);

    ssize_t nbytes_read;
    while (c != '\0') {
        // Read character
        nbytes_read = read(STDIN_FILENO, &c, 1);
        if (!nbytes_read || c == '\n')
            c = '\0';

        // Realloacte if needed
        if (str_len == buffer_size) {
            buffer_size += BLOCK_BUFFER_SIZE;
            (*pbuffer) = (char*) realloc(*pbuffer, buffer_size);
        }

        (*pbuffer)[str_len++] = c;
    }

    return str_len;
}

void main()
{
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

    ecc curve;
    ecc_init(
        curve, 
        curve_a, // a 
        curve_b, // b
        curve_p  // p
    );

    eccpt G;
    ecc_init_pt_str(
        curve, G,
        Gx,   // x
        Gy,   // y
        NULL  // z
    );

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

    mpz_t k;
    mpz_init(k);
    if (dlog(curve, k, G, kG, n, NUM_THREADS) == DLOG_SUCCESS) {
        mpz_out_str(stdout, 10, k);
        printf("\n");
    }
    else {
        printf("None\n");
    }

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