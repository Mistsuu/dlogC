#include "dlog.h"
#include "ecc.h"
#include "ecc_x.h"
#include "num.h"
#include "const.h"
#include "que2.h"
#include "mem.h"

#define mpz_size_bytes(n) mpz_sizeinbase(n, 256)

size_t dlog_calc_mem(
    mpz_t n,
    size_t* index_size_limbs, size_t* index_size_bytes,
    size_t* item_size_limbs, size_t* item_size_bytes,
    size_t* n_partitions,
    
    mpz_t upper_k,
    size_t mem_limit,
    mpz_t curve_p
)
{
    (*item_size_bytes) = mpz_size_bytes(curve_p);
    (*item_size_limbs) = mpz_size(curve_p);
    if (mem_limit && mem_limit < ((*item_size_bytes) + 1) * 4) {
        #ifdef DLOG_VERBOSE
            printf("[debug] Requested limited memory is not enough to run this algorithm. At least %ld bytes is required.\n", ((*item_size_bytes) + 1) * 4);
        #endif
        return 0;
    }

    // Number of [n | p] items we have to allocate 
    // as if we don't have memory limit.
    mpz_sqrt(n, upper_k);
    mpz_add_ui(n, n, 1);

    // Calculate number of items given the memory limit
    size_t n_limit;
    if (mem_limit) {
        for (
            size_t index_size_bytes_limit = 1, index_upper_limit = 256; 
            index_size_bytes_limit <= sizeof(size_t); 
            ++index_size_bytes_limit, index_upper_limit *= 256
        ) {
            n_limit = mem_limit / 2 / (index_size_bytes_limit + (*item_size_bytes));
            if (n_limit <= index_upper_limit - 1)
                break;
        }

        if (mpz_cmp_ui(n, n_limit - 1) > 0)
            mpz_set_ui(n, n_limit - 1);
    }

    // Set index size
    (*index_size_limbs) = mpz_size(n);
    (*index_size_bytes) = mpz_size_bytes(n);

    // Check if n fits in size_t -- argument size of malloc
    if ((*index_size_bytes) > sizeof(size_t))
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
    // size_t is argument size of malloc
    if (SIZE_MAX / ((*index_size_bytes) + (*item_size_bytes)) / 2 < n_size_t)
        return 0;

    // Set number of partitions
    mpz_t nn;
    mpz_init(nn);
    mpz_mul(nn, n, n);

    mpz_t n_partitions_q_mpz;
    mpz_t n_partitions_r_mpz;
    mpz_init(n_partitions_q_mpz);
    mpz_init(n_partitions_r_mpz);
    mpz_tdiv_qr(n_partitions_q_mpz, n_partitions_r_mpz, upper_k, nn);

    (*n_partitions) = (
         mpz_get_ui(n_partitions_q_mpz) + 
        (mpz_cmp_ui(n_partitions_r_mpz, 0) != 0)
    );
    
    mpz_clear(nn);
    mpz_clear(n_partitions_q_mpz);
    mpz_clear(n_partitions_r_mpz);
    return n_size_t;
}

int dlog_alloc_buffer(
    char** lbuffer,
    char** rbuffer,

    size_t n_size_t, 
    size_t index_size_bytes, size_t item_size_bytes
)
{
    // Checks if we're going to limit our memory?
    size_t nbytes_alloc = (n_size_t + 1) * (index_size_bytes + item_size_bytes);

    // Now we allocate
    (*lbuffer) = (char*) malloc(nbytes_alloc);
    (*rbuffer) = (char*) malloc(nbytes_alloc);
    #ifdef DLOG_VERBOSE
        printf("[debug] size buffer: %ld bytes = %f MB = %f GB\n", 
                nbytes_alloc * 2, 
                nbytes_alloc * 2 / 1024.0 / 1024.0, 
                nbytes_alloc * 2 / 1024.0 / 1024.0 / 1024.0
            );
    #endif

    // Check for valid allocation?
    if (!(*lbuffer) || !(*rbuffer)) {
        if (*lbuffer)
            free(*lbuffer);
        if (*rbuffer)
            free(*rbuffer);
        return 0;
    }

    return 1;
}

void* __thread__dlog_fill_buffer(
    void* vargs
)
{   
    
    // -------------------------------------------------------------------------------------
    //      Setup arguments.
    // -------------------------------------------------------------------------------------
    __args_thread__dlog_fill_buffer* args = (__args_thread__dlog_fill_buffer*) vargs;

    char* _buffer = args->_buffer;

    size_t n_size_t = args->n_size_t; 
    size_t index_size_bytes = args->index_size_bytes; 
    size_t index_size_limbs = args->index_size_limbs;
    size_t item_size_bytes = args->item_size_bytes; 
    size_t item_size_limbs = args->item_size_limbs;

    mp_limb_t* _0x = args->_0x; 
    mp_limb_t* _0z = args->_0z;
    mp_limb_t* _1x = args->_1x; 
    mp_limb_t* _1z = args->_1z;
    mp_limb_t* dGx = args->dGx; 
    mp_limb_t* dGz = args->dGz;

    mp_limb_t* i = args->i; 
    int is_inc_i = args->is_inc_i;

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
        // ---------------------- Write i_ -----------------------
        mpn2bytes(_buffer, index_size_bytes, i, index_size_limbs);
        _buffer += index_size_bytes;

        // ---------------------- Write _0x/_0z -----------------------
        if (ecc_xz_to_X(V0, _0x, _0z, curve_p, item_size_limbs, T))
            mpn2bytes(_buffer, item_size_bytes, V0, item_size_limbs);
        else 
            mpn2bytes(_buffer, item_size_bytes, curve_p, item_size_limbs);
        _buffer += item_size_bytes;

        // ---------------------- Update i_ -----------------------
        if (is_inc_i)
            mpn_add_1(i, i, index_size_limbs, 1);
        else
            mpn_sub_1(i, i, index_size_limbs, 1);

        // ---------------------- Update _*x,_*z -----------------------
        if (mpn_zero_p(_0z, item_size_limbs) == 0)
            ecc_xadd(
                V0, V1,
                _1x, _1z,
                dGx, dGz,
                _0x, _0z,

                curve_a,
                curve_b,
                curve_p,
                item_size_limbs,

                T
            );
        else
            ecc_xdbl(
                V0, V1,
                _1x, _1z,

                curve_a,
                curve_b,
                curve_p,
                item_size_limbs,

                T
            );
        mpn_copyd(_0x, _1x, item_size_limbs);
        mpn_copyd(_0z, _1z, item_size_limbs);
        mpn_copyd(_1x, V0, item_size_limbs);
        mpn_copyd(_1z, V1, item_size_limbs);
    }

    ecc_free_xtemp(T);
    free(V0);
    free(V1);

    return NULL;
}

void dlog_fill_buffer_l(
    char* lbuffer, 
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
    // -------------------------------------------------------------------------------------

    memset(lbuffer, 0, index_size_bytes); 
    lbuffer += index_size_bytes;
    mpn2bytes(lbuffer, item_size_bytes, mpz_limbs_read(curve->p), item_size_limbs);
    lbuffer += item_size_bytes;

    // -------------------------------------------------------------------------------------
    //      Set up L0 & L1.
    // -------------------------------------------------------------------------------------
    eccpt L0; ecc_init_pt_pt(L0, G);            // L0 <- 1*G
    eccpt L1; ecc_init_pt(L1);                  // L1 <- 2*G
    ecc_add(curve, L1, G, G);                 

    // -------------------------------------------------------------------------------------
    //      Initialize memory for arguments of __thread__dlog_fill_buffer.
    // -------------------------------------------------------------------------------------

    __args_thread__dlog_fill_buffer* thread_args = (__args_thread__dlog_fill_buffer*) malloc_exit_when_null(sizeof(__args_thread__dlog_fill_buffer) * n_threads);

    mp_limb_t* Gx = mpz_limbs_init_cpy(G->x, item_size_limbs);
    mp_limb_t* Gz = mpz_limbs_init_cpy(G->z, item_size_limbs);
    mp_limb_t* curve_a = mpz_limbs_init_cpy(curve->a, item_size_limbs);
    mp_limb_t* curve_b = mpz_limbs_init_cpy(curve->b, item_size_limbs);
    mp_limb_t* curve_p = mpz_limbs_init_cpy(curve->p, item_size_limbs);

    eccpt n_per_threadG;
    ecc_init_pt(n_per_threadG);
    ecc_mul(curve, n_per_threadG, G, n_per_thread);

    mpz_t iL;
    mpz_init_set_si(iL, 1);     // Because we have written 0 to buffer, so index of lbuffer starts with 1.

    size_t n_per_thread_size_t = n_size_t / n_threads;
    size_t n_last_thread_size_t = n_per_thread_size_t + n_size_t % n_threads;

    for (unsigned int i = 0; i < n_threads; ++i) {
        thread_args[i]._buffer = lbuffer;

        thread_args[i].n_size_t = (i != n_threads - 1 ? n_per_thread_size_t : n_last_thread_size_t);
        thread_args[i].index_size_bytes = index_size_bytes;
        thread_args[i].index_size_limbs = index_size_limbs;
        thread_args[i].item_size_bytes = item_size_bytes;
        thread_args[i].item_size_limbs = item_size_limbs;

        thread_args[i].dGx = Gx;
        thread_args[i].dGz = Gz;
        thread_args[i].curve_a = curve_a;
        thread_args[i].curve_b = curve_b;
        thread_args[i].curve_p = curve_p;

        thread_args[i].i = mpz_limbs_init_cpy(iL, index_size_limbs);
        thread_args[i].is_inc_i = 1;

        thread_args[i]._0x = mpz_limbs_init_cpy(L0->x, item_size_limbs);
        thread_args[i]._0z = mpz_limbs_init_cpy(L0->z, item_size_limbs);
        thread_args[i]._1x = mpz_limbs_init_cpy(L1->x, item_size_limbs);
        thread_args[i]._1z = mpz_limbs_init_cpy(L1->z, item_size_limbs);

        ecc_add(curve, L0, L0, n_per_threadG);      // L0 <- L0 + n_per_thread*G
        ecc_add(curve, L1, L1, n_per_threadG);      // L1 <- L1 + n_per_thread*G
        mpz_add(iL, iL, n_per_thread);

        lbuffer += n_per_thread_size_t * (item_size_bytes + index_size_bytes);
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

    ecc_free_pt(L0);
    ecc_free_pt(L1);
    ecc_free_pt(n_per_threadG);

    free(Gx);
    free(Gz);
    free(curve_a);
    free(curve_b);
    free(curve_p);

    for (unsigned int i = 0; i < n_threads; ++i) {
        free(thread_args[i].i);
        free(thread_args[i]._0x);
        free(thread_args[i]._0z);
        free(thread_args[i]._1x);
        free(thread_args[i]._1z);
    }

    free(thread_args);
    free(threads);
}

void dlog_fill_buffer_r(
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
    //     Write  [ n | (k - n*n)G ] to rbuffer.
    // -------------------------------------------------------------------------------------
    eccpt kG_nnG; ecc_init_pt_pt(kG_nnG, G);    // kG - nnG
    ecc_mul(curve, kG_nnG, kG_nnG, n);
    ecc_mul(curve, kG_nnG, kG_nnG, n);
    ecc_sub(curve, kG_nnG, kG, kG_nnG);

    // [ n | (k - n*n)G ] -> rbuffer
    mpn2bytes(rbuffer, index_size_bytes, mpz_limbs_read(n), mpz_size(n));
    rbuffer += index_size_bytes;
    mpn2bytes(rbuffer, item_size_bytes, mpz_limbs_read(kG_nnG->x), mpz_size(kG_nnG->x));
    rbuffer += item_size_bytes;

    // -------------------------------------------------------------------------------------
    //      Set up R0 & R1.
    // -------------------------------------------------------------------------------------
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

    mpz_t iR;
    mpz_init_set(iR, n);
    mpz_sub_ui(iR, iR, 1);      // Because we have written n to buffer, so index of rbuffer starts with n-1.

    size_t n_per_thread_size_t = n_size_t / n_threads;
    size_t n_last_thread_size_t = n_per_thread_size_t + n_size_t % n_threads;

    for (unsigned int i = 0; i < n_threads; ++i) {
        thread_args[i]._buffer = rbuffer;

        thread_args[i].n_size_t = (i != n_threads - 1 ? n_per_thread_size_t : n_last_thread_size_t);
        thread_args[i].index_size_bytes = index_size_bytes;
        thread_args[i].index_size_limbs = index_size_limbs;
        thread_args[i].item_size_bytes = item_size_bytes;
        thread_args[i].item_size_limbs = item_size_limbs;

        thread_args[i].dGx = nGx;
        thread_args[i].dGz = nGz;
        thread_args[i].curve_a = curve_a;
        thread_args[i].curve_b = curve_b;
        thread_args[i].curve_p = curve_p;

        thread_args[i].i = mpz_limbs_init_cpy(iR, index_size_limbs);
        thread_args[i].is_inc_i = 0;

        thread_args[i]._0x = mpz_limbs_init_cpy(R0->x, item_size_limbs);
        thread_args[i]._0z = mpz_limbs_init_cpy(R0->z, item_size_limbs);
        thread_args[i]._1x = mpz_limbs_init_cpy(R1->x, item_size_limbs);
        thread_args[i]._1z = mpz_limbs_init_cpy(R1->z, item_size_limbs);

        ecc_add(curve, R0, R0, nn_per_threadG);     // R0 <- R0 + n*n_per_thread*G
        ecc_add(curve, R1, R1, nn_per_threadG);     // R1 <- R1 + n*n_per_thread*G
        mpz_sub(iR, iR, n_per_thread);

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
    mpz_clear(iR);

    ecc_free_pt(nG);
    ecc_free_pt(kG_nnG);
    ecc_free_pt(R0);
    ecc_free_pt(R1);
    ecc_free_pt(n_per_threadG);
    ecc_free_pt(nn_per_threadG);

    free(nGx);
    free(nGz);
    free(curve_a);
    free(curve_b);
    free(curve_p);

    for (unsigned int i = 0; i < n_threads; ++i) {
        free(thread_args[i].i);
        free(thread_args[i]._0x);
        free(thread_args[i]._0z);
        free(thread_args[i]._1x);
        free(thread_args[i]._1z);
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
    char* _buffer,

    size_t n_size_t,
    size_t index_size_bytes,
    size_t item_size_bytes,
    
    unsigned int n_threads
)
{
    n_size_t += 1; // Remember, index is [0 -> n].

    size_t slot_size_bytes = index_size_bytes + item_size_bytes;
    size_t n_per_thread_size_t = n_size_t / n_threads;
    size_t n_last_thread_size_t = n_per_thread_size_t + n_size_t % n_threads;

    for (unsigned int i = 0; i < n_threads; ++i) {
        dlog_sort_one_buffer(
            _buffer, 
            i != n_threads - 1 ? n_per_thread_size_t : n_last_thread_size_t, 
            index_size_bytes, 
            item_size_bytes
        );

        _buffer += slot_size_bytes * n_per_thread_size_t;
    }
}

void* __thread__dlog_search_buffer(
    void* vargs
)
{
    // -------------------------------------------------------------------------------------
    //      Setup arguments.
    // -------------------------------------------------------------------------------------
    __args_thread__dlog_search_buffer* args = (__args_thread__dlog_search_buffer*) vargs;

    mp_limb_t* exp_l_limbs = args->exp_l_limbs;
    mp_limb_t* exp_r_limbs = args->exp_r_limbs;

    char* lbuffer = args->lbuffer;
    char* rbuffer = args->rbuffer;

    size_t n_size_t_l = args->n_size_t_l;
    size_t n_size_t_r = args->n_size_t_r;

    size_t index_size_limbs = args->index_size_limbs; 
    size_t index_size_bytes = args->index_size_bytes;
    size_t item_size_bytes = args->item_size_bytes;

    // -------------------------------------------------------------------------------------
    //      Real calculation
    // -------------------------------------------------------------------------------------

    size_t slot_size_bytes = index_size_bytes + item_size_bytes;
    char* lend = lbuffer + n_size_t_l * slot_size_bytes;
    char* rend = rbuffer + n_size_t_r * slot_size_bytes;

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
        return (void*) 0;

    mpn_zero(exp_l_limbs, index_size_limbs+1);
    mpn_zero(exp_r_limbs, index_size_limbs+1);
    mpn_set_str(exp_l_limbs, lbuffer, index_size_bytes, 256);
    mpn_set_str(exp_r_limbs, rbuffer, index_size_bytes, 256);
    return (void*) 1;
}


int dlog_search_buffer(
    mpz_t exp_l,
    mpz_t exp_r,

    char* lbuffer,
    char* rbuffer,
    
    size_t n_size_t, 
    size_t index_size_limbs, size_t index_size_bytes, 
    size_t item_size_bytes,

    unsigned int n_threads
)
{
    int is_search_found = 0;
    n_size_t += 1; // Remember, index is [0 -> n]

    // -------------------------------------------------------------------------------------
    //      Initialize memory for arguments of __thread__dlog_search_buffer.
    // -------------------------------------------------------------------------------------

    __args_thread__dlog_search_buffer* thread_args = (__args_thread__dlog_search_buffer*) malloc_exit_when_null(sizeof(__args_thread__dlog_search_buffer) * n_threads);
    
    size_t slot_size_bytes = index_size_bytes + item_size_bytes;
    size_t n_per_thread_size_t = n_size_t / n_threads;
    size_t n_last_thread_size_t = n_per_thread_size_t + n_size_t % n_threads;

    for (unsigned int i = 0; i < n_threads; ++i) {
        thread_args[i].exp_l_limbs = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * (index_size_limbs+1));
        thread_args[i].exp_r_limbs = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * (index_size_limbs+1));

        thread_args[i].index_size_bytes = index_size_bytes;
        thread_args[i].index_size_limbs = index_size_limbs;
        thread_args[i].item_size_bytes = item_size_bytes;

        thread_args[i].rbuffer = rbuffer;
        thread_args[i].n_size_t_r = i != n_threads - 1 ? n_per_thread_size_t : n_last_thread_size_t;

        rbuffer += slot_size_bytes * n_per_thread_size_t;
    }

    // -------------------------------------------------------------------------------------
    //      Generate threads.
    // -------------------------------------------------------------------------------------

    pthread_t* threads = (pthread_t*) malloc_exit_when_null(sizeof(pthread_t) * n_threads);
    int result_code;
    int is_search_found_in_thread = 0;
    for (unsigned int j = 0; j < n_threads && !is_search_found; ++j) {
        for (unsigned int i = 0; i < n_threads; ++i) {
            thread_args[i].lbuffer = lbuffer;
            thread_args[i].n_size_t_l = j != n_threads - 1 ? n_per_thread_size_t : n_last_thread_size_t;
        }

        for (unsigned int i = 0; i < n_threads; ++i) {
            result_code = pthread_create(&threads[i], NULL, __thread__dlog_search_buffer, (void*)(&thread_args[i]));
            if (result_code) {
                printf("[error] oh no! dlog_fill_buffer cannot CREATE thread!!!\n");
                exit(-1);
            }
        }

        for (unsigned int i = 0; i < n_threads; ++i) {
            result_code = pthread_join(threads[i], (void**)&is_search_found_in_thread);
            if (result_code) {
                printf("[error] oh no! dlog_fill_buffer cannot JOIN thread!!!\n");
                exit(-1);
            }

            // If result is found from search, write to exp_l :)
            if (is_search_found_in_thread && !is_search_found) {
                mpz_t tmp;
                mpz_set(exp_l, mpz_roinit_n(tmp, thread_args[i].exp_l_limbs, index_size_limbs));
                mpz_set(exp_r, mpz_roinit_n(tmp, thread_args[i].exp_r_limbs, index_size_limbs));
                is_search_found = 1;
            }
        }

        lbuffer += slot_size_bytes * n_per_thread_size_t;
    }

    // -------------------------------------------------------------------------------------
    //      Cleaning up.
    // -------------------------------------------------------------------------------------
    for (unsigned int i = 0; i < n_threads; ++i) {
        free(thread_args[i].exp_l_limbs);
        free(thread_args[i].exp_r_limbs);
    }

    free(thread_args);
    free(threads);
    return is_search_found;
}

/*
    dlog():
        ? Calculate k from G and k*G where k < upper_k.
        ! k must be init-ed before put into this function.
*/
int __dlog__(
    ecc curve, 
    mpz_t k, 
    eccpt G, eccpt kG, 
    
    char* lbuffer, 
    char* rbuffer,
    
    mpz_t n, size_t n_size_t,
    size_t index_size_limbs, size_t index_size_bytes,
    size_t item_size_limbs, size_t item_size_bytes,

    unsigned int n_threads,
    unsigned int is_update_lbuffer
)
{
    int dlog_ret_code = DLOG_NOT_FOUND_DLOG;

    mpz_t exp_l; mpz_init(exp_l);
    mpz_t exp_r; mpz_init(exp_r);
    eccpt Y; ecc_init_pt(Y);

    #ifdef DLOG_VERBOSE
        struct timeval time_start_op; 
        struct timeval time_end_op; 
        struct timeval time_elapsed_op;
    #endif

    // -------------------------------------------------------------------------------------
    //      Step 1: Filling the L & R buffers.
    // -------------------------------------------------------------------------------------

    if (is_update_lbuffer) {
        #ifdef DLOG_VERBOSE
            printf("[debug] Filling L buffer...\n");
            gettimeofday(&time_start_op, NULL);
        #endif

        dlog_fill_buffer_l(
            lbuffer, 

            curve, G, kG, 
            
            n, n_size_t, 
            index_size_bytes, index_size_limbs,
            item_size_bytes, item_size_limbs,
            
            n_threads
        );

        #ifdef DLOG_VERBOSE
            gettimeofday(&time_end_op, NULL);
            timersub(&time_end_op, &time_start_op, &time_elapsed_op);
            printf("[debug] Filling L took %ld.%06ld seconds.\n", (long int)time_elapsed_op.tv_sec, (long int)time_elapsed_op.tv_usec);
        #endif
    }

    #ifdef DLOG_VERBOSE
        printf("[debug] Filling R buffer...\n");
        gettimeofday(&time_start_op, NULL);
    #endif

    dlog_fill_buffer_r(
        rbuffer, 

        curve, G, kG, 
        
        n, n_size_t, 
        index_size_bytes, index_size_limbs,
        item_size_bytes, item_size_limbs,
        
        n_threads
    );

    #ifdef DLOG_VERBOSE
        gettimeofday(&time_end_op, NULL);
        timersub(&time_end_op, &time_start_op, &time_elapsed_op);
        printf("[debug] Filling R took %ld.%06ld seconds.\n", (long int)time_elapsed_op.tv_sec, (long int)time_elapsed_op.tv_usec);
    #endif

    // -------------------------------------------------------------------------------------
    //      Step 2: Sorting the L & R buffers.
    // -------------------------------------------------------------------------------------

    if (is_update_lbuffer) {
        #ifdef DLOG_VERBOSE
            printf("[debug] Sorting L buffer...\n");
            gettimeofday(&time_start_op, NULL);
        #endif

        dlog_sort_buffer(
            lbuffer,
        
            n_size_t, 
            index_size_bytes, 
            item_size_bytes,

            n_threads
        );

        #ifdef DLOG_VERBOSE
            gettimeofday(&time_end_op, NULL);
            timersub(&time_end_op, &time_start_op, &time_elapsed_op);
            printf("[debug] Sorting L took %ld.%06ld seconds.\n", (long int)time_elapsed_op.tv_sec, (long int)time_elapsed_op.tv_usec);
        #endif
    }

    #ifdef DLOG_VERBOSE
        printf("[debug] Sorting R buffer...\n");
        gettimeofday(&time_start_op, NULL);
    #endif

    dlog_sort_buffer(
        rbuffer,
    
        n_size_t, 
        index_size_bytes, 
        item_size_bytes,

        n_threads
    );

    #ifdef DLOG_VERBOSE
        gettimeofday(&time_end_op, NULL);
        timersub(&time_end_op, &time_start_op, &time_elapsed_op);
        printf("[debug] Sorting R took %ld.%06ld seconds.\n", (long int)time_elapsed_op.tv_sec, (long int)time_elapsed_op.tv_usec);
    #endif


    // -------------------------------------------------------------------------------------
    //      Step 3: Searching the L & R buffers.
    // -------------------------------------------------------------------------------------
    
    #ifdef DLOG_VERBOSE
        printf("[debug] Searching in L & R buffers...\n");
        gettimeofday(&time_start_op, NULL);
    #endif

    int dlog_search_status = 
        dlog_search_buffer(
            exp_l, exp_r,
        
            lbuffer,
            rbuffer,
        
            n_size_t, 
            index_size_limbs, index_size_bytes, 
            item_size_bytes,

            n_threads
        );

    #ifdef DLOG_VERBOSE
        gettimeofday(&time_end_op, NULL);
        timersub(&time_end_op, &time_start_op, &time_elapsed_op);
        printf("[debug] Searching took %ld.%06ld seconds.\n", (long int)time_elapsed_op.tv_sec, (long int)time_elapsed_op.tv_usec);
    #endif

    if (!dlog_search_status) {
        #ifdef DLOG_VERBOSE
            printf("[debug] Cannot search for equal values in L & R buffers!\n");
        #endif

        dlog_ret_code = DLOG_NOT_FOUND_DLOG;
        goto dlog_end;
    }

    // -------------------------------------------------------------------------------------
    //      Deduce k from the search.
    // -------------------------------------------------------------------------------------

    // -- Case 1: l*X = Y - r*n*X
    mpz_mul(k, exp_r, n);
    mpz_add(k, k, exp_l);
    ecc_mul(curve, Y, G, k);
    if (mpz_cmp(Y->x, kG->x) == 0 && mpz_cmp(Y->y, kG->y) == 0)
    {
        #ifdef DLOG_VERBOSE
            printf("[debug] Found k = ");
            mpz_out_str(stdout, 10, k);
            printf(".\n");
        #endif

        dlog_ret_code = DLOG_SUCCESS;
        goto dlog_end;
    }

    // -- Case 2: -l*X = Y - r*n*X
    mpz_mul(k, exp_r, n);
    mpz_sub(k, k, exp_l);
    ecc_mul(curve, Y, G, k);
    if (mpz_cmp(Y->x, kG->x) == 0 && mpz_cmp(Y->y, kG->y) == 0)
    {
        #ifdef DLOG_VERBOSE
            printf("[debug] Found k = ");
            mpz_out_str(stdout, 10, k);
            printf(".\n");
        #endif

        dlog_ret_code = DLOG_SUCCESS;
        goto dlog_end;
    }

    // -------------------------------------------------------------------------------------
    //      Cleanup.
    // -------------------------------------------------------------------------------------

dlog_end:
    mpz_clear(exp_l);
    mpz_clear(exp_r);
    ecc_free_pt(Y);

    return dlog_ret_code;
}


int dlog(
    ecc curve, 
    mpz_t k, 
    eccpt G, eccpt kG, 
    mpz_t upper_k, 

    unsigned int n_threads,
    size_t mem_limit
)
{
    #ifdef DLOG_VERBOSE
        printf("[debug] curve: \n");
        printf("[debug]    ");
        ecc_printf(curve);
        printf("\n");
        printf("[debug] G: \n");
        printf("[debug]    ");
        ecc_printf_pt(G);
        printf("\n");
        printf("[debug] kG: \n");
        printf("[debug]    ");
        ecc_printf_pt(kG);
        printf("\n");
        printf("[debug] upper_k = ");
        mpz_out_str(stdout, 10, upper_k);
        printf("\n");
        if (mem_limit)
            printf("[debug] memory limit: %ld bytes = %f MB = %f GB\n", 
                    mem_limit, 
                    mem_limit / 1024.0 / 1024.0, 
                    mem_limit / 1024.0 / 1024.0 / 1024.0
                );
        else
            printf("[debug] no memory limit is set.\n");
        printf("[debug] n_threads = %d\n", n_threads);
    #endif

    // Could have used the abs() version,
    // but this is much better.
    if (mpz_cmp_si(upper_k, 0) <= 0) {
        #ifdef DLOG_VERBOSE
            printf("[debug] Negative value of k is detected! Exiting...\n");
        #endif
        return DLOG_INVALID_UPPER_K;
    }

    // -------------------------------------------------------------------------------------
    //      Calculating the amount of memory needed for each data unit
    //      to run the algorithm.
    // -------------------------------------------------------------------------------------

    mpz_t n;
    mpz_init(n);

    size_t index_size_bytes;
    size_t item_size_bytes;
    size_t index_size_limbs;
    size_t item_size_limbs;
    size_t n_partitions;

    size_t n_size_t = dlog_calc_mem(
        n,
        &index_size_limbs, &index_size_bytes,
        &item_size_limbs, &item_size_bytes,
        &n_partitions,

        upper_k,
        mem_limit,
        curve->p
    );

    // Memory calculation failed
    if (!n_size_t) {
        #ifdef DLOG_VERBOSE
            printf("[error] Cannot allocate memory for L & R buffers!\n");
        #endif

        mpz_clear(n);
        return DLOG_CANNOT_ALLOCATE;
    }

    #ifdef DLOG_VERBOSE
        printf("[debug] index_size_bytes = %ld\n", index_size_bytes);
        printf("[debug] item_size_bytes = %ld\n", item_size_bytes);
        printf("[debug] index_size_limbs = %ld\n", index_size_limbs);
        printf("[debug] item_size_limbs = %ld\n", item_size_limbs);
        printf("[debug] n_partitions = %ld\n", n_partitions);
        printf("[debug] n_items / partition = %ld\n", n_size_t);
    #endif

    // Doing multithread this case would
    // have caused a memory error.
    if (n_size_t < n_threads) {
        #ifdef DLOG_VERBOSE
            printf("[debug] Setting n_threads = 1 because some thread will be empty...\n");
        #endif
        n_threads = 1;
    }

    // -------------------------------------------------------------------------------------
    //      Allocating the amount of memory needed to hold the buffers.
    // -------------------------------------------------------------------------------------

    char* lbuffer;
    char* rbuffer;
    if (!dlog_alloc_buffer(
        &lbuffer,
        &rbuffer,
        n_size_t, index_size_bytes, item_size_bytes
    )) 
    {
        #ifdef DLOG_VERBOSE
            printf("[error] Cannot allocate memory for L & R buffers!\n");
        #endif

        mpz_clear(n);
        return DLOG_CANNOT_ALLOCATE;
    };


    // -------------------------------------------------------------------------------------
    //      Start main operation.
    // -------------------------------------------------------------------------------------
    eccpt nnG;
    ecc_init_pt(nnG);
    ecc_mul(curve, nnG, G, n);
    ecc_mul(curve, nnG, nnG, n);

    eccpt kG_innG;
    ecc_init_pt_pt(kG_innG, kG);

    mpz_t inn;
    mpz_init(inn);
    mpz_mul(inn, n, n);

    int dlog_ret_code = DLOG_NOT_FOUND_DLOG;
    for (int n_partition = 0; n_partition < (int)n_partitions; ++n_partition) {
        #ifdef DLOG_VERBOSE
            if (n_partitions > 1) {
                printf("\n[debug] === Running partition %d/%ld === ", n_partition+1, n_partitions);
                printf("(found k value in this partition will be added with ");
                mpz_out_str(stdout, 10, inn);
                printf("*%d to get the actual k)\n", n_partition);
            }
        #endif

        dlog_ret_code = __dlog__(
            curve,
            k,
            G, kG_innG,

            lbuffer,
            rbuffer,

            n, n_size_t,
            index_size_limbs, index_size_bytes,
            item_size_limbs, item_size_bytes,

            n_threads,
            n_partition == 0
        );

        if (dlog_ret_code == DLOG_SUCCESS) {
            mpz_mul_si(inn, inn, n_partition);
            mpz_add(k, k, inn);
            break;
        }

        ecc_sub(curve, kG_innG, kG_innG, nnG);
    }

    // -------------------------------------------------------------------------------------
    //      Cleaning up.
    // -------------------------------------------------------------------------------------
    free(lbuffer);
    free(rbuffer);
    mpz_clear(n);
    mpz_clear(inn);
    ecc_free_pt(nnG);
    ecc_free_pt(kG_innG);

    return dlog_ret_code;
}