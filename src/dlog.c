#include "dlog.h"
#include "ecc.h"
#include "ecc_x.h"
#include "const.h"

#define mpz_size_bytes(n) mpz_sizeinbase(n, 256)
            
size_t dlog_alloc_buffer(
    char** buffer,
    mpz_t n, size_t index_size_bytes, size_t item_size_bytes
)
{
    // Check if n fits in size_t -- argument size of malloc
    if (mpz_size_bytes(n) > sizeof(size_t))
        return 0;
    printf("[debug] mpz_size_bytes(n) = %ld\n", mpz_size_bytes(n));

    // Convert n from mpz_t to size_t
    size_t           nitems_alloc_size   = mpz_size(n);
    const mp_limb_t* nitems_alloc_limbs  = mpz_limbs_read(n);
    size_t           n_size_t            = 0;
    for (int i = nitems_alloc_size-1; i >= 0; --i) {
        n_size_t <<= mp_bits_per_limb;
        n_size_t ^=  nitems_alloc_limbs[i];
    }

    // Checks if n_size_t items fits in size_t
    if (SIZE_MAX / (index_size_bytes + item_size_bytes) / 2 < n_size_t)
        return 0;

    // Now we allocate!
    size_t nbytes_alloc = n_size_t * (index_size_bytes + item_size_bytes) * 2 + 1;
    printf("[debug] size buffer: %ld bytes = %f MB = %f GB\n", 
                nbytes_alloc, 
                nbytes_alloc / 1024.0 / 1024.0, 
                nbytes_alloc / 1024.0 / 1024.0 / 1024.0
            );
    (*buffer) = (char*) malloc(nbytes_alloc);
    return n_size_t;
}

void __thread__dlog_fill_buffer(
    char* buffer,

    mpz_t n, size_t n_size_t, 
    size_t index_size_bytes, size_t index_size_limbs,
    size_t item_size_bytes, size_t item_size_limbs,

    mp_limb_t* Lx, mp_limb_t* Lz,   // must be modified in-use
    mp_limb_t* Rx, mp_limb_t* Rz,   // must be modified in-use
    mp_limb_t* Px, mp_limb_t* Pz,

    mp_limb_t* i_item_l, mp_limb_t* i_item_r,

    mp_limb_t* curve_a,
    mp_limb_t* curve_b,
    mp_limb_t* curve_p
)
{

}

void dlog_fill_buffer(
    char* buffer, 
    ecc curve, eccpt G, eccpt kG, 

    mpz_t n, size_t n_size_t, 
    size_t index_size_bytes, size_t index_size_limbs,
    size_t item_size_bytes, size_t item_size_limbs,
    
    unsigned int n_threads
)
{
    
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
    mpz_t n;
    mpz_init(n);
    mpz_sqrt(n, upper_k);
    mpz_add_ui(n, n, 1);
    
    size_t index_size_bytes = mpz_size_bytes(n);
    size_t item_size_bytes  = mpz_size_bytes(curve->p);
    size_t index_size_limbs = mpz_size(n);
    size_t item_size_limbs  = mpz_size(curve->p);

    char* buffer;
    size_t n_size_t = dlog_alloc_buffer(
        &buffer,
        n, index_size_bytes, item_size_bytes
    );

    // Allocation failed
    if (!n_size_t) {
        mpz_clear(n);
        printf("[debug] cannot allocate memory!\n");
        return DLOG_CANNOT_ALLOCATE;
    }
    printf("[debug] index_size_bytes = %ld\n", index_size_bytes);
    printf("[debug] item_size_bytes = %ld\n", item_size_bytes);

    dlog_fill_buffer(
        buffer, 
        curve, G, kG, 
        
        n, n_size_t, 
        index_size_bytes, index_size_limbs,
        item_size_bytes, item_size_limbs,
        
        n_threads
    );
    // dlog_sort_buffer(
    //     buffer, 
    //     n_size_t, index_size_bytes, item_size_bytes,
    //     n_threads
    // );

    // mpz_t iL, iR;
    // dlog_search_buffer(
    //     iL, iR,
    //     buffer, 
    //     n_size_t, index_size_bytes, item_size_bytes,
    //     n_threads
    // );

    mpz_clear(n);
    free(buffer);
    return DLOG_SUCCESS;
}