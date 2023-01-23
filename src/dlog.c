#include "dlog.h"
#include "ecc.h"
#include "ecc_x.h"
#include "const.h"

#define mpz_size_bytes(n) (                             \
            ((mpz_size(n) * mp_bits_per_limb) >> 3) +    \
            ((mpz_size(n) * mp_bits_per_limb) &  7 != 0) \
        )                                               \

mp_limb_t* dlog_allocate(mpz_t nitems_alloc, size_t n_size_limbs, size_t p_size_limbs)
{
    // Check if nitems_alloc fits in size_t -- argument size of malloc
    if (mpz_size_bytes(nitems_alloc) > sizeof(size_t))
        return NULL;

    // Convert nitems_alloc from mpz_t to size_t
    size_t           nitems_alloc_size   = mpz_size(nitems_alloc);
    const mp_limb_t* nitems_alloc_limbs  = mpz_limbs_read(nitems_alloc);
    size_t           size_t_nitems_alloc = 0;
    for (int i = nitems_alloc_size-1; i >= 0; --i) {
        size_t_nitems_alloc <<= mp_bits_per_limb;
        size_t_nitems_alloc ^=  nitems_alloc_limbs[i];
    }

    // Checks if sizeof(mp_limb_t) * size_t_nitems_alloc fits in size_t
    if (SIZE_MAX / sizeof(mp_limb_t) / (n_size_limbs + p_size_limbs) < size_t_nitems_alloc)
        return NULL;

    // Now we allocate!
    size_t nbytes_alloc = sizeof(mp_limb_t) * size_t_nitems_alloc * (n_size_limbs + p_size_limbs);
    printf("[i] Debug size buffer: %ld bytes = %f MB = %f GB\n", 
                nbytes_alloc, 
                nbytes_alloc / 1024.0 / 1024.0, 
                nbytes_alloc / 1024.0 / 1024.0 / 1024.0
            );
    return (mp_limb_t*) malloc(sizeof(mp_limb_t) * size_t_nitems_alloc * (n_size_limbs + p_size_limbs));
}

/*
    dlog():
        ? Calculate k from G and k*G where k < n.
        ! k must be init-ed before put into this function.
*/
int dlog(ecc curve, mpz_t k, eccpt G, eccpt kG, mpz_t n)
{
    assert(mpz_cmp_si(n, 4) > 0);
    size_t n_size_limbs = mpz_size(n);
    size_t p_size_limbs = mpz_size(curve->p);
    
    // Number of [n | p] items we have to allocate.
    mpz_t nitems_alloc;
    mpz_init(nitems_alloc);
    mpz_sqrt(nitems_alloc, n);
    mpz_add_ui(nitems_alloc, nitems_alloc, 1);

    mp_limb_t* buffer = dlog_allocate(nitems_alloc, n_size_limbs, p_size_limbs);
    if (!buffer) {
        mpz_clear(nitems_alloc);
        return DLOG_CANNOT_ALLOCATE;
    }

    mpz_clear(nitems_alloc);
    free(buffer);
    return DLOG_SUCCESS;
}