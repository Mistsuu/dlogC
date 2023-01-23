#include "dlog.h"
#include "ecc.h"
#include "ecc_x.h"
#include "num.h"
#include "const.h"

#define mpz_size_bytes(n) mpz_sizeinbase(n, 256)
            
size_t dlog_init_buffer(
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

    // Now we allocate (and set 0 the whole area :>)
    size_t nbytes_alloc = n_size_t * (index_size_bytes + item_size_bytes) * 2 + 1;
    printf("[debug] size buffer: %ld bytes = %f MB = %f GB\n", 
                nbytes_alloc, 
                nbytes_alloc / 1024.0 / 1024.0, 
                nbytes_alloc / 1024.0 / 1024.0 / 1024.0
            );
    (*buffer) = (char*) malloc(nbytes_alloc);
    memset(*buffer, 0, nbytes_alloc);
    return n_size_t;
}

void __thread__dlog_fill_buffer(
    char* buffer,

    size_t n_size_t, 
    size_t index_size_bytes, size_t index_size_limbs,
    size_t item_size_bytes, size_t item_size_limbs,

    mp_limb_t* L0x, mp_limb_t* L0z,   // must be modified in-use
    mp_limb_t* L1x, mp_limb_t* L1z,   // must be modified in-use
    mp_limb_t* R0x, mp_limb_t* R0z,   // must be modified in-use
    mp_limb_t* R1x, mp_limb_t* R1z,   // must be modified in-use
    mp_limb_t* Px, mp_limb_t* Pz,

    mp_limb_t* i_item_l, mp_limb_t* i_item_r,

    mp_limb_t* curve_a,
    mp_limb_t* curve_b,
    mp_limb_t* curve_p
)
{   
    ecc_xtemp T;
    ecc_init_xtemp(T, item_size_limbs);

    mp_limb_t* V0 = (mp_limb_t*) malloc(sizeof(mp_limb_t) * 2 * item_size_limbs);
    mp_limb_t* V1 = (mp_limb_t*) malloc(sizeof(mp_limb_t) * item_size_limbs);
    mp_limb_t* V2 = (mp_limb_t*) malloc(sizeof(mp_limb_t) * mpn_sec_invert_itch(item_size_limbs));

    for (size_t _ = 0; _ < n_size_t; ++_) {
        // ---------------------- Write iL -----------------------
        mpn_get_str(buffer, 256, i_item_l, index_size_limbs);
        buffer += index_size_bytes;

        // ---------------------- Write L1x/L1z -----------------------
        mpn_copyd(V0, L1z, item_size_limbs);
        if (mpn_sec_invert(V1, V0, curve_p, item_size_limbs, 2 * item_size_bytes * 8, V2)) {            // V1 <- 1/L1z mod p
            mpn_mul_n(V0, L1x, V1, item_size_limbs);                                                    // V0 <- L1x/L1z
            mpn_tdiv_qr(V1, V0, 0, V0, 2*item_size_limbs, curve_p, item_size_limbs);                    // V0 <- L1x/L1z mod p
            mpn_get_str(buffer, 256, V0, item_size_limbs);                                              // Write V0 to buffer.
        }
        else {
            mpn_get_str(buffer, 256, curve_p, item_size_limbs);
        }
        buffer += item_size_limbs;

        // ---------------------- Write iR -----------------------
        mpn_get_str(buffer, 256, i_item_r, index_size_limbs);
        buffer += index_size_bytes;

        // ---------------------- Write R1x/R1z -----------------------
        mpn_copyd(V0, R1z, item_size_limbs);
        if (mpn_sec_invert(V1, V0, curve_p, item_size_limbs, 2 * item_size_bytes * 8, V2)) {             // V1 <- 1/R1z mod p
            mpn_mul_n(V0, R1x, V1, item_size_limbs);                                                     // V0 <- R1x/R1z
            mpn_tdiv_qr(V1, V0, 0, V0, 2*item_size_limbs, curve_p, item_size_limbs);                     // V0 <- R1x/R1z mod p
            mpn_get_str(buffer, 256, V0, item_size_limbs);                                               // Write V0 to buffer.
        }
        else {
            mpn_get_str(buffer, 256, curve_p, item_size_limbs);
        }
        buffer += item_size_limbs;

        // ---------------------- Update iL -----------------------
        mpn_add_1(i_item_l, i_item_l, item_size_limbs, 1);

        // ---------------------- Update iR -----------------------
        mpn_sub_1(i_item_r, i_item_r, item_size_limbs, 1);

        // ---------------------- Update L*x,L*z -----------------------
        ecc_xadd(
            V0, V1,
            L1x, L1z,
            Px, Pz,
            L0x, L0z,

            curve_a,
            curve_b,
            curve_p,
            item_size_limbs,

            T
        );
        mpn_copyd(L0x, L1x, item_size_limbs);
        mpn_copyd(L0z, L1z, item_size_limbs);
        mpn_copyd(L1x, V0, item_size_limbs);
        mpn_copyd(L1z, V1, item_size_limbs);

        // ---------------------- Update R*x,R*z -----------------------
        ecc_xadd(
            V0, V1,
            R1x, R1z,
            Px, Pz,
            R0x, R0z,

            curve_a,
            curve_b,
            curve_p,
            item_size_limbs,

            T
        );
        mpn_copyd(R0x, R1x, item_size_limbs);
        mpn_copyd(R0z, R1z, item_size_limbs);
        mpn_copyd(R1x, V0, item_size_limbs);
        mpn_copyd(R1z, V1, item_size_limbs);
    }

    ecc_free_xtemp(T);
    free(V0);
    free(V1);
    free(V2);
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
    mpz_t _0;
    mpz_init_set_ui(_0, 0);

    mpz_t n_per_thread;
    mpz_init_set(n_per_thread, n);
    mpz_tdiv_r_ui(n_per_thread, n_per_thread, (unsigned long)n_threads);

    eccpt L0; ecc_init_pt_pt(L0, G);    // L0 <- G
    eccpt L1; ecc_init_pt_pt(L1, G);    
    ecc_add(curve, L1, L1, L1);         // L1 <- 2*G

    eccpt _nnG; ecc_init_pt_pt(_nnG, G);
    ecc_mul(curve, _nnG, _nnG, n);
    ecc_mul(curve, _nnG, _nnG, n);
    mpz_sub(_nnG->y, curve->p, _nnG->y);

    eccpt R0; ecc_init_pt_pt(R0, kG);   // R0 <- (k-n*n)*G
    ecc_add(curve, R0, R0, _nnG);
    eccpt R1; ecc_init_pt_pt(R1, R0);    
    ecc_add(curve, R1, R1, G);          // R1 <- (k-n*n+1)*G

    mp_limb_t* iL = mpz_limbs_init_cpy(_0,           index_size_limbs);
    mp_limb_t* iR = mpz_limbs_init_cpy(n_per_thread, index_size_limbs);
    mp_limb_t* Gx = mpz_limbs_init_cpy(G->x, item_size_limbs);
    mp_limb_t* Gz = mpz_limbs_init_cpy(G->z, item_size_limbs);
    mp_limb_t* L0x = mpz_limbs_init_cpy(L0->x, item_size_limbs);
    mp_limb_t* L0z = mpz_limbs_init_cpy(L0->z, item_size_limbs);
    mp_limb_t* L1x = mpz_limbs_init_cpy(L1->x, item_size_limbs);
    mp_limb_t* L1z = mpz_limbs_init_cpy(L1->z, item_size_limbs);
    mp_limb_t* R0x = mpz_limbs_init_cpy(R0->x, item_size_limbs);
    mp_limb_t* R0z = mpz_limbs_init_cpy(R0->z, item_size_limbs);
    mp_limb_t* R1x = mpz_limbs_init_cpy(R1->x, item_size_limbs);
    mp_limb_t* R1z = mpz_limbs_init_cpy(R1->z, item_size_limbs);
    mp_limb_t* curve_a = mpz_limbs_init_cpy(curve->a, item_size_limbs);
    mp_limb_t* curve_b = mpz_limbs_init_cpy(curve->b, item_size_limbs);
    mp_limb_t* curve_p = mpz_limbs_init_cpy(curve->p, item_size_limbs);

    size_t n_per_thread_size_t = n_size_t / n_threads;
    __thread__dlog_fill_buffer(
        buffer,

        n_per_thread_size_t, 
        index_size_bytes, index_size_limbs,
        item_size_bytes, item_size_limbs,

        L0x, L0z,
        L1x, L1z,
        R0x, R0z,
        R1x, R1z,
        Gx, Gz,

        iL, iR,

        curve_a,
        curve_b,
        curve_p
    );

    mpz_clear(_0);
    mpz_clear(n_per_thread);

    ecc_free_pt(L0);
    ecc_free_pt(L1);
    ecc_free_pt(R0);
    ecc_free_pt(R1);

    free(Gx);
    free(Gz);
    free(L0x);
    free(L0z);
    free(L1x);
    free(L1z);
    free(R0x);
    free(R0z);
    free(R1x);
    free(R1z);
    free(curve_a);
    free(curve_b);
    free(curve_p);
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
    size_t n_size_t = dlog_init_buffer(
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