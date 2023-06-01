#include "rand.h"
#include "mem.h"

#define mpz_size_bytes(n) mpz_sizeinbase(n, 256)

void mpz_dev_urandomm(mpz_t rop, const mpz_t n)
{
    FILE *fptr = fopen("/dev/urandom", "r");
    if (!fptr) {
        printf("[error] Error! Cannot read from random source!\n");
        exit(-1);
    }

    // Read random bytes
    // to mp_limb_t.
    size_t item_size_limbs = mpz_size(n);
    mp_limb_t* tmp_limbs; 
    tmp_limbs = (mp_limb_t*) malloc_exit_when_null(item_size_limbs * sizeof(mp_limb_t));
    if (fread(tmp_limbs, sizeof(mp_limb_t), item_size_limbs, fptr) != item_size_limbs) {
        printf("[error] Error! Memory source is insufficient!\n");
        exit(-1);
    }

    // Convert to mpz_t
    mpz_t tmp;
    mpz_set(rop, mpz_roinit_n(tmp, tmp_limbs, item_size_limbs));
    mpz_mod(rop, rop, n);

    // Free memory.
    free(tmp_limbs);
    fclose(fptr);
}