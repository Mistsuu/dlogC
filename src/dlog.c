#include "dlog.h"
#include "ecc.h"
#include "ecc_x.h"
#include "const.h"

#define mpz_size_bytes(n) mpz_sizeinbase(n, 256)
            
char* dlog_alloc_buffer(mpz_t nitems, size_t index_size, size_t item_size)
{
    // Check if nitems fits in size_t -- argument size of malloc
    if (mpz_size_bytes(nitems) > sizeof(size_t))
        return NULL;

    printf("[debug] mpz_size_bytes(nitems) = %ld\n", mpz_size_bytes(nitems));

    // Convert nitems from mpz_t to size_t
    size_t           nitems_alloc_size   = mpz_size(nitems);
    const mp_limb_t* nitems_alloc_limbs  = mpz_limbs_read(nitems);
    size_t           size_t_nitems_alloc = 0;
    for (int i = nitems_alloc_size-1; i >= 0; --i) {
        size_t_nitems_alloc <<= mp_bits_per_limb;
        size_t_nitems_alloc ^=  nitems_alloc_limbs[i];
    }

    // Checks if size_t_nitems_alloc items fits in size_t
    if (SIZE_MAX / (index_size + item_size) / 2 < size_t_nitems_alloc)
        return NULL;

    // Now we allocate!
    size_t nbytes_alloc = size_t_nitems_alloc * (index_size + item_size) * 2 + 1;
    printf("[debug] size buffer: %ld bytes = %f MB = %f GB\n", 
                nbytes_alloc, 
                nbytes_alloc / 1024.0 / 1024.0, 
                nbytes_alloc / 1024.0 / 1024.0 / 1024.0
            );
    return (char*) malloc(nbytes_alloc);
}

/*
    dlog():
        ? Calculate k from G and k*G where k < upper_k.
        ! k must be init-ed before put into this function.
*/
int dlog(ecc curve, mpz_t k, eccpt G, eccpt kG, mpz_t upper_k, unsigned int n_threads)
{
    assert(mpz_cmp_si(upper_k, 4) > 0);
    
    // Number of [n | p] items we have to allocate.
    mpz_t nitems;
    mpz_init(nitems);
    mpz_sqrt(nitems, upper_k);
    mpz_add_ui(nitems, nitems, 1);
    
    size_t index_size = mpz_size_bytes(nitems);
    size_t item_size  = mpz_size_bytes(curve->p);
    char* buffer = dlog_alloc_buffer(nitems, index_size, item_size);
    if (!buffer) {
        mpz_clear(nitems);
        printf("[debug] cannot allocate memory!\n");
        return DLOG_CANNOT_ALLOCATE;
    }

    printf("[debug] index_size = %ld\n", index_size);
    printf("[debug] item_size = %ld\n", item_size);

    dlog_fill_buffer(
        buffer, 
        curve, G, kG, 
        nitems, index_size, item_size,
        n_threads
    );
    // dlog_sort_buffer(
    //     buffer, 
    //     nitems, index_size, item_size,
    //     n_threads
    // );

    // mpz_t iL, iR;
    // dlog_search_buffer(
    //     iL, iR,
    //     buffer, 
    //     nitems, index_size, item_size,
    //     n_threads
    // );

    mpz_clear(nitems);
    free(buffer);
    return DLOG_SUCCESS;
}