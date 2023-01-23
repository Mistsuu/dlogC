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
    mpz_set_str(mpz_a, "18446744073709551618", 10);
    mpz_set_str(mpz_m, "18446744073709551619", 10);

    // Get input space
    mp_limb_t* mpn_a = NULL;
    mp_limb_t* mpn_m = NULL;
    mp_size_t mpn_input_size = get_mp_input_limps(&mpn_a, &mpn_m, mpz_a, mpz_m);

    printf("[i] A mem read:\n");
    for (int i = 0; i < mpn_input_size; ++i) {
        printf(" L %d: 0x%016lx\n", i, mpn_a[i]);
    }

    printf("[i] M mem read:\n");
    for (int i = 0; i < mpn_input_size; ++i) {
        printf(" L %d: 0x%016lx\n", i, mpn_m[i]);
    }

    mpz_clear(mpz_a);
    mpz_clear(mpz_m);
    free(mpn_a);
    free(mpn_m);
}