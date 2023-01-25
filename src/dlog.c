#include "dlog.h"
#include "ecc.h"
#include "ecc_x.h"
#include "num.h"
#include "const.h"
#include "que2.h"
#include "mem.h"

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
    (*lbuffer) = (char*) malloc_exit_when_null(nbytes_alloc);
    (*rbuffer) = (char*) malloc_exit_when_null(nbytes_alloc);
    return n_size_t;
}

void* __thread__dlog_fill_buffer(
    void* vargs
)
{   
    
    // -------------------------------------------------------------------------------------
    //      Setup arguments.
    // -------------------------------------------------------------------------------------
    __args_thread__dlog_fill_buffer* args = (__args_thread__dlog_fill_buffer*) vargs;

    char* lbuffer = args->lbuffer;
    char* rbuffer = args->rbuffer;

    size_t n_size_t = args->n_size_t; 
    size_t index_size_bytes = args->index_size_bytes; 
    size_t index_size_limbs = args->index_size_limbs;
    size_t item_size_bytes = args->item_size_bytes; 
    size_t item_size_limbs = args->item_size_limbs;

    mp_limb_t* L0x = args->L0x; 
    mp_limb_t* L0z = args->L0z;
    mp_limb_t* L1x = args->L1x; 
    mp_limb_t* L1z = args->L1z;
    mp_limb_t* R0x = args->R0x; 
    mp_limb_t* R0z = args->R0z;
    mp_limb_t* R1x = args->R1x; 
    mp_limb_t* R1z = args->R1z;
    mp_limb_t* Gx = args->Gx; 
    mp_limb_t* Gz = args->Gz;
    mp_limb_t* nGx = args->nGx; 
    mp_limb_t* nGz = args->nGz;

    mp_limb_t* i_l = args->i_l; 
    mp_limb_t* i_r = args->i_r;

    mp_limb_t* curve_a = args->curve_a;
    mp_limb_t* curve_b = args->curve_b;
    mp_limb_t* curve_p = args->curve_p;

    // -------------------------------------------------------------------------------------
    //      Real calculation
    // -------------------------------------------------------------------------------------

    ecc_xtemp T;
    ecc_init_xtemp(T, item_size_limbs);

    mp_limb_t* V0 = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * item_size_limbs);
    mp_limb_t* V1 = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * item_size_limbs);

    for (size_t _ = 0; _ < n_size_t; ++_) {
        // ---------------------- Write iL -----------------------
        mpn2bytes(lbuffer, index_size_bytes, i_l, index_size_limbs);
        lbuffer += index_size_bytes;

        // ---------------------- Write L0x/L0z -----------------------
        if (ecc_xz_to_X(V0, L0x, L0z, curve_p, item_size_limbs, T))
            mpn2bytes(lbuffer, item_size_bytes, V0, item_size_limbs);
        else 
            mpn2bytes(lbuffer, item_size_bytes, curve_p, item_size_limbs);
        lbuffer += item_size_bytes;

        // ---------------------- Write iR -----------------------
        mpn2bytes(rbuffer, index_size_bytes, i_r, index_size_limbs);
        rbuffer += index_size_bytes;

        // ---------------------- Write R0x/R0z -----------------------
        if (ecc_xz_to_X(V0, R0x, R0z, curve_p, item_size_limbs, T))
            mpn2bytes(rbuffer, item_size_bytes, V0, item_size_limbs);
        else 
            mpn2bytes(rbuffer, item_size_bytes, curve_p, item_size_limbs);
        rbuffer += item_size_bytes;

        // ---------------------- Update iL -----------------------
        mpn_add_1(i_l, i_l, index_size_limbs, 1);

        // ---------------------- Update iR -----------------------
        mpn_sub_1(i_r, i_r, index_size_limbs, 1);

        // ---------------------- Update L*x,L*z -----------------------
        ecc_xadd(
            V0, V1,
            L1x, L1z,
            Gx, Gz,
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
            nGx, nGz,
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

    __args_thread__dlog_fill_buffer* thread_args = (__args_thread__dlog_fill_buffer*) malloc_exit_when_null(sizeof(__args_thread__dlog_fill_buffer) * n_threads);

    mp_limb_t* Gx = mpz_limbs_init_cpy(G->x, item_size_limbs);
    mp_limb_t* Gz = mpz_limbs_init_cpy(G->z, item_size_limbs);
    mp_limb_t* nGx = mpz_limbs_init_cpy(nG->x, item_size_limbs);
    mp_limb_t* nGz = mpz_limbs_init_cpy(nG->z, item_size_limbs);
    mp_limb_t* curve_a = mpz_limbs_init_cpy(curve->a, item_size_limbs);
    mp_limb_t* curve_b = mpz_limbs_init_cpy(curve->b, item_size_limbs);
    mp_limb_t* curve_p = mpz_limbs_init_cpy(curve->p, item_size_limbs);

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

    size_t n_per_thread_size_t = n_size_t / n_threads;
    size_t n_last_thread_size_t = n_per_thread_size_t + n_size_t % n_threads;

    for (unsigned int i = 0; i < n_threads; ++i) {
        thread_args[i].lbuffer = lbuffer;
        thread_args[i].rbuffer = rbuffer;
        thread_args[i].n_size_t = (i != n_threads - 1 ? n_per_thread_size_t : n_last_thread_size_t);
        thread_args[i].index_size_bytes = index_size_bytes;
        thread_args[i].index_size_limbs = index_size_limbs;
        thread_args[i].item_size_bytes = item_size_bytes;
        thread_args[i].item_size_limbs = item_size_limbs;

        thread_args[i].Gx = Gx;
        thread_args[i].nGx = nGx;
        thread_args[i].Gz = Gz;
        thread_args[i].nGz = nGz;
        thread_args[i].curve_a = curve_a;
        thread_args[i].curve_b = curve_b;
        thread_args[i].curve_p = curve_p;

        thread_args[i].i_l = mpz_limbs_init_cpy(iL, index_size_limbs);
        thread_args[i].i_r = mpz_limbs_init_cpy(iR, index_size_limbs);
        thread_args[i].L0x = mpz_limbs_init_cpy(L0->x, item_size_limbs);
        thread_args[i].L0z = mpz_limbs_init_cpy(L0->z, item_size_limbs);
        thread_args[i].L1x = mpz_limbs_init_cpy(L1->x, item_size_limbs);
        thread_args[i].L1z = mpz_limbs_init_cpy(L1->z, item_size_limbs);
        thread_args[i].R0x = mpz_limbs_init_cpy(R0->x, item_size_limbs);
        thread_args[i].R0z = mpz_limbs_init_cpy(R0->z, item_size_limbs);
        thread_args[i].R1x = mpz_limbs_init_cpy(R1->x, item_size_limbs);
        thread_args[i].R1z = mpz_limbs_init_cpy(R1->z, item_size_limbs);

        ecc_add(curve, L0, L0, n_per_threadG);      // L0 <- L0 + n_per_thread*G
        ecc_add(curve, L1, L1, n_per_threadG);      // L1 <- L1 + n_per_thread*G
        ecc_add(curve, R0, R0, nn_per_threadG);     // R0 <- R0 + n*n_per_thread*G
        ecc_add(curve, R1, R1, nn_per_threadG);     // R1 <- R1 + n*n_per_thread*G
        mpz_add(iL, iL, n_per_thread);
        mpz_sub(iR, iR, n_per_thread);

        lbuffer += n_per_thread_size_t * (item_size_bytes + index_size_bytes);
        rbuffer += n_per_thread_size_t * (item_size_bytes + index_size_bytes);
    }

    // -------------------------------------------------------------------------------------
    //      Generate threads.
    // -------------------------------------------------------------------------------------

    pthread_t* threads = (pthread_t*) malloc_exit_when_null(sizeof(pthread_t) * n_threads);
    int result_code;
    for (unsigned int i = 0; i < n_threads; ++i) {
        result_code = pthread_create(&threads[i], NULL, __thread__dlog_fill_buffer, (void*)(&thread_args[i]));
        if (result_code) {
            printf("[error] oh no! dlog_fill_buffer cannot CREATE thread!!!\n");
            exit(-1);
        }
    }

    for (unsigned int i = 0; i < n_threads; ++i) {
        result_code = pthread_join(threads[i], NULL);
        if (result_code) {
            printf("[error] oh no! dlog_fill_buffer cannot JOIN thread!!!\n");
            exit(-1);
        }
    }

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
        free(thread_args[i].i_l);
        free(thread_args[i].i_r);
        free(thread_args[i].L0x);
        free(thread_args[i].L0z);
        free(thread_args[i].L1x);
        free(thread_args[i].L1z);
        free(thread_args[i].R0x);
        free(thread_args[i].R0z);
        free(thread_args[i].R1x);
        free(thread_args[i].R1z);
    }

    free(thread_args);
    free(threads);
}

void dlog_sort_one_buffer(
    char* buffer,

    size_t n_size_t,
    size_t index_size_bytes,
    size_t item_size_bytes
)
{
    // Create a temporary buffer to have a place holder for swapping :)
    size_t slot_size_bytes = index_size_bytes + item_size_bytes;
    char* temp_buf = (char*) malloc_exit_when_null(slot_size_bytes);

    // Queue for holding lo, hi values.
    que2 sort_queue;
    que2_init(sort_queue);

    signed long LO;
    signed long HI;
    signed long PI;

    char* pivot;
    char* buffer_i = buffer;
    char* buffer_j = buffer;

    que2_push(sort_queue, 0, n_size_t-1);
    while (que2_pop(sort_queue, &LO, &HI)) {
        pivot    = buffer + HI * slot_size_bytes;
        buffer_i = buffer + LO * slot_size_bytes;
        buffer_j = buffer + LO * slot_size_bytes;
        PI       = LO;
        while (buffer_j < pivot) {
            // Swap i, j views if element is < pivot
            if (memcmp(
                    buffer_j + index_size_bytes, 
                    pivot + index_size_bytes, 
                    item_size_bytes
                ) < 0) {
                    // Swap i view with j view.
                    memcpy(temp_buf, buffer_i, slot_size_bytes);
                    memcpy(buffer_i, buffer_j, slot_size_bytes);
                    memcpy(buffer_j, temp_buf, slot_size_bytes);

                    // Update i view.
                    buffer_i += slot_size_bytes;
                    PI       += 1;
                }

            // Update j view.
            buffer_j += slot_size_bytes;
        }

        // Swap i view with pivot :)
        memcpy(temp_buf, buffer_i, slot_size_bytes);
        memcpy(buffer_i, pivot, slot_size_bytes);
        memcpy(pivot, temp_buf, slot_size_bytes);

        // Append lookup indices to quicksort queue.
        if (LO < PI-1)
            que2_push(sort_queue, LO, PI-1);
        if (PI+1 < HI)
            que2_push(sort_queue, PI+1, HI);
    }

    que2_free(sort_queue);
    free(temp_buf);
}

void dlog_sort_buffer(
    char* lbuffer,
    char* rbuffer,

    size_t n_size_t,
    size_t index_size_bytes,
    size_t item_size_bytes
)
{
    n_size_t += 1; // Remember, index is [0 -> n].
    dlog_sort_one_buffer(lbuffer, n_size_t, index_size_bytes, item_size_bytes);
    dlog_sort_one_buffer(rbuffer, n_size_t, index_size_bytes, item_size_bytes);
}


int dlog_search_buffer(
    mpz_t exp_l,
    mpz_t exp_r,

    char* lbuffer,
    char* rbuffer,
    
    size_t n_size_t, 
    size_t index_size_limbs, size_t index_size_bytes, 
    size_t item_size_bytes
)
{
    size_t slot_size_bytes = index_size_bytes + item_size_bytes;
    char* lend = lbuffer + n_size_t * slot_size_bytes;
    char* rend = rbuffer + n_size_t * slot_size_bytes;

    int cmp_status;
    while (lbuffer != lend && rbuffer != rend) {
        cmp_status = memcmp(
                        lbuffer + index_size_bytes, 
                        rbuffer + index_size_bytes, 
                        item_size_bytes
                     );

        if (cmp_status < 0)
            lbuffer += slot_size_bytes;
        else if (cmp_status > 0)
            rbuffer += slot_size_bytes;
        else
            break;
    }

    if (lbuffer == lend || rbuffer == rend)
        return 0;
    
    mp_limb_t* exp_l_limbs = (mp_limb_t*) malloc(sizeof(mp_limb_t) * (index_size_limbs+1));
    mp_limb_t* exp_r_limbs = (mp_limb_t*) malloc(sizeof(mp_limb_t) * (index_size_limbs+1));
    mpn_zero(exp_l_limbs, index_size_limbs+1);
    mpn_zero(exp_r_limbs, index_size_limbs+1);
    printf("%ld\n", mpn_set_str(exp_l_limbs, lbuffer, index_size_bytes, 256));
    printf("%ld\n", mpn_set_str(exp_r_limbs, rbuffer, index_size_bytes, 256));

    mpz_t tmp;
    mpz_set(exp_l, mpz_roinit_n(tmp, exp_l_limbs, index_size_limbs));
    mpz_set(exp_r, mpz_roinit_n(tmp, exp_r_limbs, index_size_limbs));

    free(exp_l_limbs);
    free(exp_r_limbs);
    return 1;
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

    printf("[debug] Filling lbuffer - rbuffer...\n");
    dlog_fill_buffer(
        lbuffer, 
        rbuffer, 

        curve, G, kG, 
        
        n, n_size_t, 
        index_size_bytes, index_size_limbs,
        item_size_bytes, item_size_limbs,
        
        n_threads
    );

    printf("[debug] Sorting lbuffer - rbuffer...\n");
    dlog_sort_buffer(
        lbuffer,
        rbuffer,
    
        n_size_t, 
        index_size_bytes, 
        item_size_bytes
    );

    printf("[debug] Searching lbuffer - rbuffer...\n");
    mpz_t exp_l; mpz_init(exp_l);
    mpz_t exp_r; mpz_init(exp_r);
    if (!dlog_search_buffer(
        exp_l, exp_r,
    
        lbuffer,
        rbuffer,
    
        n_size_t, 
        index_size_limbs, index_size_bytes, 
        item_size_bytes
    )) 
    {
        mpz_clear(n);
        mpz_clear(exp_l);
        mpz_clear(exp_r);

        free(lbuffer);
        free(rbuffer);

        return DLOG_NOT_FOUND_DLOG;
    }

    // FILE* pfile;
    // pfile = fopen("miscellaneous/outputl", "w");
    // fwrite(lbuffer, 1, (n_size_t + 1) * (index_size_bytes + item_size_bytes), pfile);
    // fclose(pfile);
    // pfile = fopen("miscellaneous/outputr", "w");
    // fwrite(rbuffer, 1, (n_size_t + 1) * (index_size_bytes + item_size_bytes), pfile);
    // fclose(pfile);

    eccpt Y;
    ecc_init_pt(Y);

    // -- Case 1: l*X = Y - r*n*X
    mpz_mul(k, exp_r, n);
    mpz_add(k, k, exp_l);
    mpz_mod(k, k, curve->p);
    ecc_mul(curve, Y, G, k);
    if (mpz_cmp(Y->x, kG->x) == 0 && mpz_cmp(Y->y, kG->y) == 0)
    {
        ecc_free_pt(Y);

        mpz_clear(n);
        mpz_clear(exp_l);
        mpz_clear(exp_r);

        free(lbuffer);
        free(rbuffer);

        return DLOG_SUCCESS;
    }

    // -- Case 2: -l*X = Y - r*n*X
    mpz_mul(k, exp_r, n);
    mpz_sub(k, k, exp_l);
    mpz_mod(k, k, curve->p);
    ecc_mul(curve, Y, G, k);
    if (mpz_cmp(Y->x, kG->x) == 0 && mpz_cmp(Y->y, kG->y) == 0)
    {
        ecc_free_pt(Y);

        mpz_clear(n);
        mpz_clear(exp_l);
        mpz_clear(exp_r);

        free(lbuffer);
        free(rbuffer);

        return DLOG_SUCCESS;
    }

    ecc_free_pt(Y);

    mpz_clear(n);
    mpz_clear(exp_l);
    mpz_clear(exp_r);

    free(lbuffer);
    free(rbuffer);
    return DLOG_NOT_FOUND_DLOG;
}