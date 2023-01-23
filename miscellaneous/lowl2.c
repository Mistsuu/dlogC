#include "gmp.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

mp_size_t get_mp_input_limps(mp_limb_t** mp_limb_out1, mp_limb_t** mp_limb_out2, mpz_t mpz_in1, mpz_t mpz_in2)
{
    mp_size_t mp_biggest_size_in = mpz_in1->_mp_size;
    if (mpz_in2->_mp_size > mpz_in1->_mp_size)
        mp_biggest_size_in = mpz_in2->_mp_size;

    (*mp_limb_out1) = malloc(sizeof(mp_limb_t) * mp_biggest_size_in);
    (*mp_limb_out2) = malloc(sizeof(mp_limb_t) * mp_biggest_size_in);

    memset(*mp_limb_out1, 0, sizeof(mp_limb_t) * mp_biggest_size_in);
    memset(*mp_limb_out2, 0, sizeof(mp_limb_t) * mp_biggest_size_in);

    memcpy(*mp_limb_out1, mpz_in1->_mp_d, sizeof(mp_limb_t) * mpz_in1->_mp_size);
    memcpy(*mp_limb_out2, mpz_in2->_mp_d, sizeof(mp_limb_t) * mpz_in2->_mp_size);

    return mp_biggest_size_in;
}

int main(int argc, char * argv[]) {

    mpz_t mpz_a;
    mpz_t mpz_m;
    mpz_init(mpz_a);
    mpz_init(mpz_m);
    mpz_set_str(mpz_a, "111111111111111111111111111111111111111", 10);
    mpz_set_str(mpz_m, "211111111111111111111111111111111111111", 10);

    // Get input space
    mp_limb_t* mpn_a = NULL;
    mp_limb_t* mpn_m = NULL;
    mp_size_t mpn_input_size = get_mp_input_limps(&mpn_a, &mpn_m, mpz_a, mpz_m);
    
    // Get output space
    mp_limb_t* mpn_r = (mp_limb_t*) malloc(sizeof(mp_limb_t) * mpn_input_size);
    memset(mpn_r, 0, sizeof(mp_limb_t) * mpn_input_size);

    // Calculate
    mp_limb_t borrow = mpn_sub_n(mpn_r, mpn_a, mpn_m, mpn_input_size);
    printf("[i] Borrowing: %ld\n", borrow);

    mp_limb_t carry = mpn_add_n(mpn_r, mpn_r, mpn_m, mpn_input_size);
    printf("[i] Carrying: %ld\n", carry);

    // Print number
    mpz_t mpz_r;
    mpz_r->_mp_alloc = mpn_input_size;
    mpz_r->_mp_size  = mpn_input_size;
    mpz_r->_mp_d     = mpn_r;
    mpz_out_str(stdout, 10, mpz_r);
    printf("\n");

    printf("[i] Mem read:\n");
    for (int i = 0; i < mpn_input_size; ++i) {
        printf(" L %d: 0x%016lx\n", i, mpn_r[i]);
    }

    mpz_clear(mpz_a);
    mpz_clear(mpz_m);
    free(mpn_a);
    free(mpn_m);
    free(mpn_r);
}