#include "dlog.h"
#include "ecc.h"
#include "ecc_x.h"
#include "num.h"
#include "const.h"

#define mpz_size_bytes(n) mpz_sizeinbase(n, 256)
            
size_t dlog_init_buffer(
    char** lbuffer,
    char** rbuffer,
    mpz_t n, size_t index_size_bytes, size_t item_size_bytes
)
{
    // Check if n fits in size_t -- argument size of malloc
    if (mpz_size_bytes(n) > sizeof(size_t))
        return 0;

    // Convert n from mpz_t to size_t
    size_t           nitems_alloc_size   = mpz_size(n);
    const mp_limb_t* nitems_alloc_limbs  = mpz_limbs_read(n);
    size_t           n_size_t            = 0;
    for (int i = nitems_alloc_size-1; i >= 0; --i) {
        n_size_t <<= mp_bits_per_limb;
        n_size_t ^=  nitems_alloc_limbs[i];
    }

    // Checks if size of n_size_t items fits in size_t
    if (SIZE_MAX / (index_size_bytes + item_size_bytes) / 2 < n_size_t)
        return 0;

    // Now we allocate
    size_t nbytes_alloc = (n_size_t + 1) * (index_size_bytes + item_size_bytes);
    printf("[debug] size buffer: %ld bytes = %f MB = %f GB\n", 
                nbytes_alloc * 2, 
                nbytes_alloc * 2 / 1024.0 / 1024.0, 
                nbytes_alloc * 2 / 1024.0 / 1024.0 / 1024.0
            );
    (*lbuffer) = (char*) malloc(nbytes_alloc);
    (*rbuffer) = (char*) malloc(nbytes_alloc);
    return n_size_t;
}

void __thread__dlog_fill_buffer(
    __args_thread__dlog_fill_buffer* args
)
{   
    ecc_xtemp T;
    ecc_init_xtemp(T, args->item_size_limbs);

    mp_limb_t* V0 = (mp_limb_t*) malloc(sizeof(mp_limb_t) * args->item_size_limbs);
    mp_limb_t* V1 = (mp_limb_t*) malloc(sizeof(mp_limb_t) * args->item_size_limbs);

    for (size_t _ = 0; _ < args->n_size_t; ++_) {
        // ---------------------- Write iL -----------------------
        mpn2bytes(args->lbuffer, args->index_size_bytes, args->i_l, args->index_size_limbs);
        args->lbuffer += args->index_size_bytes;

        // ---------------------- Write L0x/L0z -----------------------
        if (ecc_xz_to_X(V0, args->L0x, args->L0z, args->curve_p, args->item_size_limbs, T))
            mpn2bytes(args->lbuffer, args->item_size_bytes, V0, args->item_size_limbs);
        else 
            mpn2bytes(args->lbuffer, args->item_size_bytes, args->curve_p, args->item_size_limbs);
        args->lbuffer += args->item_size_bytes;

        // ---------------------- Write iR -----------------------
        mpn2bytes(args->rbuffer, args->index_size_bytes, args->i_r, args->index_size_limbs);
        args->rbuffer += args->index_size_bytes;

        // ---------------------- Write R0x/R0z -----------------------
        if (ecc_xz_to_X(V0, args->R0x, args->R0z, args->curve_p, args->item_size_limbs, T))
            mpn2bytes(args->rbuffer, args->item_size_bytes, V0, args->item_size_limbs);
        else 
            mpn2bytes(args->rbuffer, args->item_size_bytes, args->curve_p, args->item_size_limbs);
        args->rbuffer += args->item_size_bytes;

        // ---------------------- Update iL -----------------------
        mpn_add_1(args->i_l, args->i_l, args->index_size_limbs, 1);

        // ---------------------- Update iR -----------------------
        mpn_sub_1(args->i_r, args->i_r, args->index_size_limbs, 1);

        // ---------------------- Update L*x,L*z -----------------------
        ecc_xadd(
            V0, V1,
            args->L1x, args->L1z,
            args->Gx, args->Gz,
            args->L0x, args->L0z,

            args->curve_a,
            args->curve_b,
            args->curve_p,
            args->item_size_limbs,

            T
        );
        mpn_copyd(args->L0x, args->L1x, args->item_size_limbs);
        mpn_copyd(args->L0z, args->L1z, args->item_size_limbs);
        mpn_copyd(args->L1x, V0, args->item_size_limbs);
        mpn_copyd(args->L1z, V1, args->item_size_limbs);

        // ---------------------- Update R*x,R*z -----------------------
        ecc_xadd(
            V0, V1,
            args->R1x, args->R1z,
            args->nGx, args->nGz,
            args->R0x, args->R0z,

            args->curve_a,
            args->curve_b,
            args->curve_p,
            args->item_size_limbs,

            T
        );
        mpn_copyd(args->R0x, args->R1x, args->item_size_limbs);
        mpn_copyd(args->R0z, args->R1z, args->item_size_limbs);
        mpn_copyd(args->R1x, V0, args->item_size_limbs);
        mpn_copyd(args->R1z, V1, args->item_size_limbs);
    }

    ecc_free_xtemp(T);
    free(V0);
    free(V1);
}

void dlog_fill_buffer(
    char* lbuffer, 
    char* rbuffer, 
    ecc curve, eccpt G, eccpt kG, 

    mpz_t n, size_t n_size_t, 
    size_t index_size_bytes, size_t index_size_limbs,
    size_t item_size_bytes, size_t item_size_limbs,
    
    unsigned int n_threads
)
{
    // -------------------------------------------------------------------------------------
    //      Set up number of items per thread.
    // -------------------------------------------------------------------------------------
    mpz_t n_per_thread;
    mpz_t n_last_thread;
    mpz_init_set(n_per_thread, n);
    mpz_init(n_last_thread);
    mpz_tdiv_qr_ui(n_per_thread, n_last_thread, n_per_thread, (unsigned long)n_threads);
    mpz_add(n_last_thread, n_last_thread, n_per_thread);

    // -------------------------------------------------------------------------------------
    //      Write [ 0 |     0*G    ] to lbuffer, 
    //            [ n | (k - n*n)G ] to rbuffer.
    // -------------------------------------------------------------------------------------
    eccpt kG_nnG; ecc_init_pt_pt(kG_nnG, G);    // kG - nnG
    ecc_mul(curve, kG_nnG, kG_nnG, n);
    ecc_mul(curve, kG_nnG, kG_nnG, n);
    ecc_sub(curve, kG_nnG, kG, kG_nnG);

    // [ 0 |     0*G    ] -> lbuffer
    memset(lbuffer, 0, index_size_bytes); 
    lbuffer += index_size_bytes;
    mpn2bytes(lbuffer, item_size_bytes, mpz_limbs_read(curve->p), item_size_limbs);
    lbuffer += item_size_bytes;

    // [ n | (k - n*n)G ] -> rbuffer
    mpn2bytes(rbuffer, index_size_bytes, mpz_limbs_read(n), index_size_limbs);
    rbuffer += index_size_bytes;
    mpn2bytes(rbuffer, item_size_bytes, mpz_limbs_read(kG_nnG->x), item_size_limbs);
    rbuffer += item_size_bytes;

    // -------------------------------------------------------------------------------------
    //      Set up L0 & L1.
    // -------------------------------------------------------------------------------------
    eccpt L0; ecc_init_pt_pt(L0, G);            // L0 <- 1*G
    eccpt L1; ecc_init_pt(L1);                  // L1 <- 2*G
    ecc_add(curve, L1, G, G);                 

    eccpt nG; ecc_init_pt_pt(nG, G);
    ecc_mul(curve, nG, nG, n);
    eccpt R0; ecc_init_pt_pt(R0, kG_nnG);       // R0 <- (k - n*n + 1*n)*G
    ecc_add(curve, R0, R0, nG);
    eccpt R1; ecc_init_pt_pt(R1, R0);    
    ecc_add(curve, R1, R1, nG);                 // R1 <- (k - n*n + 2*n)*G

    // -------------------------------------------------------------------------------------
    //      Initialize memory for arguments of __thread__dlog_fill_buffer.
    // -------------------------------------------------------------------------------------

    __args_thread__dlog_fill_buffer* thread_args = (__args_thread__dlog_fill_buffer*) malloc(sizeof(__args_thread__dlog_fill_buffer) * n_threads);

    mp_limb_t* Gx = mpz_limbs_init_cpy(G->x, item_size_limbs);
    mp_limb_t* Gz = mpz_limbs_init_cpy(G->z, item_size_limbs);
    mp_limb_t* nGx = mpz_limbs_init_cpy(nG->x, item_size_limbs);
    mp_limb_t* nGz = mpz_limbs_init_cpy(nG->z, item_size_limbs);
    mp_limb_t* curve_a = mpz_limbs_init_cpy(curve->a, item_size_limbs);
    mp_limb_t* curve_b = mpz_limbs_init_cpy(curve->b, item_size_limbs);
    mp_limb_t* curve_p = mpz_limbs_init_cpy(curve->p, item_size_limbs);

    mp_limb_t** i_l = (mp_limb_t**) malloc(sizeof(mp_limb_t*) * n_threads);
    mp_limb_t** i_r = (mp_limb_t**) malloc(sizeof(mp_limb_t*) * n_threads);
    mp_limb_t** L0x = (mp_limb_t**) malloc(sizeof(mp_limb_t*) * n_threads);
    mp_limb_t** L0z = (mp_limb_t**) malloc(sizeof(mp_limb_t*) * n_threads);
    mp_limb_t** L1x = (mp_limb_t**) malloc(sizeof(mp_limb_t*) * n_threads);
    mp_limb_t** L1z = (mp_limb_t**) malloc(sizeof(mp_limb_t*) * n_threads);
    mp_limb_t** R0x = (mp_limb_t**) malloc(sizeof(mp_limb_t*) * n_threads);
    mp_limb_t** R0z = (mp_limb_t**) malloc(sizeof(mp_limb_t*) * n_threads);
    mp_limb_t** R1x = (mp_limb_t**) malloc(sizeof(mp_limb_t*) * n_threads);
    mp_limb_t** R1z = (mp_limb_t**) malloc(sizeof(mp_limb_t*) * n_threads);

    eccpt n_per_threadG;
    eccpt nn_per_threadG;
    ecc_init_pt(n_per_threadG);
    ecc_init_pt(nn_per_threadG);
    ecc_mul(curve, n_per_threadG, G, n_per_thread);
    ecc_mul(curve, nn_per_threadG, n_per_threadG, n);

    mpz_t iL;
    mpz_t iR;
    mpz_init_set_si(iL, 1);     // Because we have written 0 to buffer, so index of lbuffer starts with 1.
    mpz_init_set(iR, n);
    mpz_sub_ui(iR, iR, 1);      // Because we have written n to buffer, so index of rbuffer starts with n-1.

    for (unsigned int i = 0; i < n_threads; ++i) {
        i_l[i] = mpz_limbs_init_cpy(iL, index_size_limbs);
        i_r[i] = mpz_limbs_init_cpy(iR, index_size_limbs);
        L0x[i] = mpz_limbs_init_cpy(L0->x, item_size_limbs);
        L0z[i] = mpz_limbs_init_cpy(L0->z, item_size_limbs);
        L1x[i] = mpz_limbs_init_cpy(L1->x, item_size_limbs);
        L1z[i] = mpz_limbs_init_cpy(L1->z, item_size_limbs);
        R0x[i] = mpz_limbs_init_cpy(R0->x, item_size_limbs);
        R0z[i] = mpz_limbs_init_cpy(R0->z, item_size_limbs);
        R1x[i] = mpz_limbs_init_cpy(R1->x, item_size_limbs);
        R1z[i] = mpz_limbs_init_cpy(R1->z, item_size_limbs);

        ecc_add(curve, L0, L0, n_per_threadG);      // L0 <- L0 + n_per_thread*G
        ecc_add(curve, L1, L1, n_per_threadG);      // L1 <- L1 + n_per_thread*G
        ecc_add(curve, R0, R0, nn_per_threadG);     // R0 <- R0 + n*n_per_thread*G
        ecc_add(curve, R1, R1, nn_per_threadG);     // R1 <- R1 + n*n_per_thread*G
        mpz_add(iL, iL, n_per_thread);
        mpz_sub(iR, iR, n_per_thread);
    }

    // -------------------------------------------------------------------------------------
    //      Generate threads.
    // -------------------------------------------------------------------------------------

    size_t n_per_thread_size_t = n_size_t / n_threads;
    size_t n_last_thread_size_t = n_per_thread_size_t + n_size_t % n_threads;
    // __thread__dlog_fill_buffer(
    //     lbuffer,
    //     rbuffer,

    //     n_per_thread_size_t, 
    //     index_size_bytes, index_size_limbs,
    //     item_size_bytes, item_size_limbs,

    //     L0x, L0z,
    //     L1x, L1z,
    //     R0x, R0z,
    //     R1x, R1z,
    //     Gx, Gz,
    //     nGx, nGz,

    //     i_l, i_r,

    //     curve_a,
    //     curve_b,
    //     curve_p
    // );


    // -------------------------------------------------------------------------------------
    //      Cleaning up.
    // -------------------------------------------------------------------------------------

    mpz_clear(n_per_thread);
    mpz_clear(n_last_thread);
    mpz_clear(iL);
    mpz_clear(iR);

    ecc_free_pt(nG);
    ecc_free_pt(kG_nnG);
    ecc_free_pt(L0);
    ecc_free_pt(L1);
    ecc_free_pt(R0);
    ecc_free_pt(R1);
    ecc_free_pt(n_per_threadG);
    ecc_free_pt(nn_per_threadG);

    free(Gx);
    free(Gz);
    free(nGx);
    free(nGz);
    free(curve_a);
    free(curve_b);
    free(curve_p);

    for (unsigned int i = 0; i < n_threads; ++i) {
        free(i_l[i]);
        free(i_r[i]);
        free(L0x[i]);
        free(L0z[i]);
        free(L1x[i]);
        free(L1z[i]);
        free(R0x[i]);
        free(R0z[i]);
        free(R1x[i]);
        free(R1z[i]);
    }

    free(i_l);
    free(i_r);
    free(L0x);
    free(L0z);
    free(L1x);
    free(L1z);
    free(R0x);
    free(R0z);
    free(R1x);
    free(R1z);
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
    
    size_t index_size_bytes = mpz_size_bytes(n);
    size_t item_size_bytes  = mpz_size_bytes(curve->p);
    size_t index_size_limbs = mpz_size(n);
    size_t item_size_limbs  = mpz_size(curve->p);

    char* lbuffer;
    char* rbuffer;
    size_t n_size_t = dlog_init_buffer(
        &lbuffer,
        &rbuffer,
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
    printf("[debug] index_size_limbs = %ld\n", index_size_limbs);
    printf("[debug] item_size_limbs = %ld\n", item_size_limbs);

    dlog_fill_buffer(
        lbuffer, 
        rbuffer, 

        curve, G, kG, 
        
        n, n_size_t, 
        index_size_bytes, index_size_limbs,
        item_size_bytes, item_size_limbs,
        
        n_threads
    );

    // dlog_sort_buffer(
    //     lbuffer,
    //     rbuffer,
    //
    //     n_size_t, 
    //     index_size_bytes, 
    //     item_size_bytes,
    //
    //     n_threads
    // );

    // mpz_t iL, iR;
    // dlog_search_buffer(
    //     iL, iR,
    //
    //     lbuffer,
    //     rbuffer,
    //
    //     n_size_t, 
    //     index_size_bytes, 
    //     item_size_bytes,
    //
    //     n_threads
    // );

    mpz_clear(n);
    free(lbuffer);
    free(rbuffer);
    
    return DLOG_SUCCESS;
}