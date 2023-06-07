#include "dlog.h"
#include "ecc.h"
#include "ecc_proj.h"
#include "ex_mpn.h"
#include "ex_mpz.h"
#include "const.h"
#include "mem.h"

int dlog_validate_input(
    mpz_t p,
    mpz_t G, mpz_t kG,
    mpz_t G_mult_order,

    unsigned int n_threads,
    unsigned int n_cache_items,
    unsigned int n_rand_items
)
{
    // -------------------------------------------------------------------------------------
    //      Must have more than 0 threads.
    //      Cache must not have 0 items.
    //      Random indices must not have 0 items.
    // -------------------------------------------------------------------------------------
    if (n_threads == 0) {
        #ifdef DLOG_VERBOSE
            printf("[debug] 0 threads is also not supported. Exiting...\n");
        #endif
        return DLOG_BAD_CONFIG;
    }

    if (n_cache_items == 0) {
        #ifdef DLOG_VERBOSE
            printf("[debug] You must have >= 1 item in the cache. Exiting...\n");
        #endif
        return DLOG_BAD_CONFIG;
    }

    if (n_rand_items < 2) {
        #ifdef DLOG_VERBOSE
            printf("[debug] You must have >= 2 random elements. Exiting...\n");
        #endif
        return DLOG_BAD_CONFIG;
    }

    // -------------------------------------------------------------------------------------
    //      Reduce G, kG mod p
    // -------------------------------------------------------------------------------------
    mpz_mod(G, G, p);
    mpz_mod(kG, kG, p);

    // -------------------------------------------------------------------------------------
    //      Is p prime?
    // -------------------------------------------------------------------------------------
    if (!(mpz_probab_prime_p(p, BASESIZE_PRIME_CHECKER))) {
        #ifdef DLOG_VERBOSE
            printf("[debug] p isn't prime. Exiting...\n");
        #endif
        return DLOG_P_IS_NOT_PRIME;
    }

    // -------------------------------------------------------------------------------------
    //      Is G's order positive and prime?
    // -------------------------------------------------------------------------------------
    if (mpz_cmp_si(G_mult_order, 0) <= 0) {
        #ifdef DLOG_VERBOSE
            printf("[debug] G_mult_order is negative. Exiting...\n");
        #endif
        return DLOG_FAULTY_POINT_ORDER;
    }

    if (!mpz_probab_prime_p(G_mult_order, BASESIZE_PRIME_CHECKER)) {
        #ifdef DLOG_VERBOSE
            printf("[debug] G_mult_order isn't prime. Exiting...\n");
        #endif
        return DLOG_FAULTY_POINT_ORDER;
    }

    // -------------------------------------------------------------------------------------
    //      Check G^Gorder = 1.
    //      Not checking kG^Gorder = 1 here since
    //      it's not a pre-requisite for this algorithm
    //      to run.
    // -------------------------------------------------------------------------------------
    mpz_t O;
    mpz_init(O);
    mpz_powm(O, G, G_mult_order, p);

    if (mpz_cmp_ui(O, 1) != 0) {
        #ifdef DLOG_VERBOSE
            printf("[debug] G's order isn't G_mult_order. Exiting...\n");
        #endif
        mpz_clear(O);
        return DLOG_FAULTY_POINT_ORDER;
    }

    mpz_clear(O);
    return DLOG_MOVE_TO_NEXT_STEP;
}

int dlog_fast_solve_if_possible(
    mpz_t p,
    mpz_t k, 
    mpz_t G, mpz_t kG,
    mpz_t G_mult_order
)
{
    // -------------------------------------------------------------------------------------
    //      Case 0: kG = 1
    // -------------------------------------------------------------------------------------
    if (mpz_cmp_ui(kG, 1) == 0) {
            #ifdef DLOG_VERBOSE
            printf("[debug] kG == 1. Automatically set k = 0.\n");
        #endif
        mpz_set_ui(k, 0);
        return DLOG_SUCCESS;
    }

    // -------------------------------------------------------------------------------------
    //      Case 1: G == 1 and kG != 1
    // -------------------------------------------------------------------------------------
    if (mpz_cmp_ui(G, 1) == 0) {
        #ifdef DLOG_VERBOSE
            printf("[debug] G == 1 and kG != 1. No solution exists.\n");
        #endif
        return DLOG_NOT_FOUND_DLOG;
    }

    // -------------------------------------------------------------------------------------
    //      Case 2: kG.order() > G.order()
    // -------------------------------------------------------------------------------------
    mpz_t O;
    mpz_init(O);
    mpz_powm(O, kG, G_mult_order, p);
    if (mpz_cmp_ui(O, 1) != 0) {
        #ifdef DLOG_VERBOSE
            printf("[debug] G^n == 1 but kG^n != 1. No solution exists.\n");
        #endif
        mpz_clear(O);
        return DLOG_NOT_FOUND_DLOG;
    }

    // -------------------------------------------------------------------------------------
    //      Case 3: I don't know any function
    //              that can help us to determine that
    //              kG is from G^k. So this should be
    //              left here blank.
    // -------------------------------------------------------------------------------------
    mpz_clear(O);
    return DLOG_MOVE_TO_NEXT_STEP;
}

void dlog_init_dlog_obj(
    dlog_obj obj,

    mpz_t p,
    mpz_t G, mpz_t kG,
    mpz_t G_mult_order,

    unsigned int n_threads,
    unsigned int n_cache_items,
    unsigned int n_rand_items
)
{
    obj->n_threads = n_threads;
    obj->n_cache_items = n_cache_items;
    obj->n_rand_items = n_rand_items;
    obj->item_size_limbs  = mpz_size(p);
    obj->index_size_limbs = mpz_size(G_mult_order);

    obj->field_p = mpn_init_zero(obj->item_size_limbs);
    obj->field_P = mpn_init_zero(obj->item_size_limbs);
    obj->G_order = mpn_init_zero(obj->item_size_limbs);

    obj->thread_tortoise_items = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_threads);
    obj->thread_tortoise_ts_indices = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_threads);
    obj->thread_hare_items_caches = (mp_limb_t***)malloc_exit_when_null(sizeof(mp_limb_t**) * n_threads);
    obj->thread_hare_ts_index_caches = (mp_limb_t***)malloc_exit_when_null(sizeof(mp_limb_t**) * n_threads);

    for (unsigned int ithread = 0; ithread < n_threads; ++ithread) {
        obj->thread_tortoise_items[ithread] = mpn_init_zero(obj->item_size_limbs);
        obj->thread_tortoise_ts_indices[ithread] = mpn_init_zero(obj->index_size_limbs * 2);
        
        obj->thread_hare_items_caches[ithread] = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_cache_items);
        obj->thread_hare_ts_index_caches[ithread] = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_cache_items);
        for (unsigned int jthread = 0; jthread < n_cache_items; ++jthread) {
            obj->thread_hare_items_caches[ithread][jthread] = mpn_init_zero(obj->item_size_limbs);
            obj->thread_hare_ts_index_caches[ithread][jthread] = mpn_init_zero(obj->index_size_limbs * 2);
        }
    }

    obj->thread_read_counters = (unsigned long***)malloc_exit_when_null(sizeof(unsigned long**) * n_threads);
    obj->thread_write_counters = (unsigned long**)malloc_exit_when_null(sizeof(unsigned long*) * n_threads);

    for (unsigned int ithread = 0; ithread < n_threads; ++ithread) {
        obj->thread_read_counters[ithread] = (unsigned long**)malloc_exit_when_null(sizeof(unsigned long*) * n_threads);
        obj->thread_write_counters[ithread] = (unsigned long*)malloc_exit_when_null(sizeof(unsigned long) * n_cache_items);
        for (unsigned int jthread = 0; jthread < n_threads; ++jthread) {
            obj->thread_read_counters[ithread][jthread] = (unsigned long*)malloc_exit_when_null(sizeof(unsigned long) * n_cache_items);
        }
    }

    obj->random_tG_add_skG = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_rand_items);
    obj->random_ts         = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_rand_items);
    for (unsigned int irand = 0; irand < n_rand_items; ++irand) {
        obj->random_tG_add_skG[irand] = mpn_init_zero(obj->item_size_limbs);
        obj->random_ts[irand] = mpn_init_zero(obj->index_size_limbs * 2);
    }

    obj->thread_result_tortoise_ts_indices = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_threads);
    obj->thread_result_hare_ts_indices     = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_threads);
    for (unsigned int ithread = 0; ithread < n_threads; ++ithread) {
        obj->thread_result_tortoise_ts_indices[ithread] = mpn_init_zero(obj->index_size_limbs * 2);
        obj->thread_result_hare_ts_indices[ithread]     = mpn_init_zero(obj->index_size_limbs * 2);
    }

    obj->founds = (int*)malloc_exit_when_null(sizeof(int) * n_threads);

    #ifdef DLOG_VERBOSE
        obj->thread_cache_hit_counters              = (mpz_t*)malloc_exit_when_null(sizeof(mpz_t) * n_threads);
        obj->thread_cache_miss_counters             = (mpz_t*)malloc_exit_when_null(sizeof(mpz_t) * n_threads);
        obj->thread_cache_possible_misread_counters = (mpz_t*)malloc_exit_when_null(sizeof(mpz_t) * n_threads);

        for (unsigned int i = 0; i < n_threads; ++i) {
            mpz_init(obj->thread_cache_hit_counters[i]);
            mpz_init(obj->thread_cache_miss_counters[i]);
            mpz_init(obj->thread_cache_possible_misread_counters[i]);
        }
    #endif
}

void dlog_fill_dlog_obj(
    dlog_obj obj,

    mpz_t p,
    mpz_t G, mpz_t kG,
    mpz_t G_mult_order,

    unsigned int n_threads,
    unsigned int n_cache_items,
    unsigned int n_rand_items
)
{
    mpz_t mpz_R;
    mpz_init_set_ui(mpz_R, 1);
    mpz_mul_2exp(mpz_R, mpz_R, obj->item_size_limbs * mp_bits_per_limb);

    // -------------------------------------------------------------------------------------
    //      Initialize p and p^-1 mod R.
    //      (very important for multiplication optimization...)
    // -------------------------------------------------------------------------------------
    mpz_t mpz_P;
    mpz_init(mpz_P);
    mpz_invert(mpz_P, p, mpz_R);
    mpz_sub(mpz_P, mpz_R, mpz_P);

    mpn_cpyz(obj->field_p, p,     obj->item_size_limbs);
    mpn_cpyz(obj->field_P, mpz_P, obj->item_size_limbs);
    mpn_cpyz(obj->G_order, G_mult_order, obj->item_size_limbs);

    // -------------------------------------------------------------------------------------
    //      Initialize fixed random elements
    // -------------------------------------------------------------------------------------
    mpz_t tG;
    mpz_t skG;
    mpz_t tG_add_skG;
    mpz_t t;
    mpz_t s;
    mpz_init(t);
    mpz_init(s);
    mpz_init(tG);
    mpz_init(skG);
    mpz_init(tG_add_skG);
    
    for (unsigned int irand = 0; irand < n_rand_items; ++irand) {
        do {
            mpz_dev_urandomm(t, G_mult_order);
            mpz_dev_urandomm(s, G_mult_order);

            mpz_powm(tG, G, t, p);
            mpz_powm(skG, kG, s, p);
            mpz_mul(tG_add_skG, tG, skG);
            mpz_mod(tG_add_skG, tG_add_skG, p);
        } while (mpz_cmp_ui(tG_add_skG, 1) == 0);

        // Convert to Montgomery form.
        mpz_mul(tG_add_skG, tG_add_skG, mpz_R);

        // Fill t,s values.
        mpn_cpyz( obj->random_ts[irand],                        t, obj->index_size_limbs);
        mpn_cpyz(&obj->random_ts[irand][obj->index_size_limbs], s, obj->index_size_limbs);
        mpn_cpyz( obj->random_tG_add_skG[irand], tG_add_skG, obj->item_size_limbs);
    }

    // -------------------------------------------------------------------------------------
    //      Initialize cache values
    // -------------------------------------------------------------------------------------
    for (unsigned int ithread = 0; ithread < n_threads; ++ithread) {
        do {
            mpz_dev_urandomm(t, G_mult_order);
            mpz_dev_urandomm(s, G_mult_order);

            mpz_powm(tG, G, t, p);
            mpz_powm(skG, kG, s, p);
            mpz_mul(tG_add_skG, tG, skG);
            mpz_mod(tG_add_skG, tG_add_skG, p);
        } while (mpz_cmp_ui(tG_add_skG, 1) == 0);

        // Convert to Montgomery form.
        mpz_mul(tG_add_skG, tG_add_skG, mpz_R);

        mpn_cpyz( obj->thread_tortoise_items[ithread], tG_add_skG, obj->item_size_limbs);
        mpn_cpyz( obj->thread_tortoise_ts_indices[ithread],                        t, obj->index_size_limbs);
        mpn_cpyz(&obj->thread_tortoise_ts_indices[ithread][obj->index_size_limbs], s, obj->index_size_limbs);

        for (unsigned int icache = 0; icache < n_cache_items; ++icache) {
            mpn_cpyz( obj->thread_hare_items_caches[ithread][icache], tG_add_skG, obj->item_size_limbs);
            mpn_cpyz( obj->thread_hare_ts_index_caches[ithread][icache],                        t, obj->index_size_limbs);
            mpn_cpyz(&obj->thread_hare_ts_index_caches[ithread][icache][obj->index_size_limbs], s, obj->index_size_limbs);
        }
    }

    // -------------------------------------------------------------------------------------
    //      Initialize read/write counters
    // -------------------------------------------------------------------------------------
    for (unsigned int ithread = 0; ithread < n_threads; ++ithread) {
        memset(obj->thread_write_counters[ithread], 0, n_cache_items * sizeof(unsigned long));
        for (unsigned int jthread = 0; jthread < n_threads; ++jthread)
            memset(obj->thread_read_counters[ithread][jthread], 0, n_cache_items * sizeof(unsigned long));
    }

    // -------------------------------------------------------------------------------------
    //      Initialize overall results
    // -------------------------------------------------------------------------------------
    memset(obj->founds, 0, sizeof(int) * n_threads);
    obj->overall_found = 0;

    // -------------------------------------------------------------------------------------
    //      Initialize some profiling variables.
    // -------------------------------------------------------------------------------------
    #ifdef DLOG_VERBOSE
        for (unsigned int ithread = 0; ithread < n_threads; ++ithread) {
            mpz_set_ui(obj->thread_cache_hit_counters[ithread], 0);
            mpz_set_ui(obj->thread_cache_miss_counters[ithread], 0);
            mpz_set_ui(obj->thread_cache_possible_misread_counters[ithread], 0);
        }
    #endif

    // -------------------------------------------------------------------------------------
    //      Free stuffs
    // -------------------------------------------------------------------------------------
    mpz_clear(mpz_P);
    mpz_clear(mpz_R);
    mpz_clear(t);
    mpz_clear(s);
    mpz_clear(tG);
    mpz_clear(skG);
    mpz_clear(tG_add_skG);
}

void dlog_free_dlog_obj(
    dlog_obj obj
)
{
    free(obj->field_p);
    free(obj->field_P);
    free(obj->G_order);

    for (unsigned int ithread = 0; ithread < obj->n_threads; ++ithread) {
        free(obj->thread_tortoise_items[ithread]);
        free(obj->thread_tortoise_ts_indices[ithread]);

        for (unsigned int icache = 0; icache < obj->n_cache_items; ++icache) {
            free(obj->thread_hare_items_caches[ithread][icache]);
            free(obj->thread_hare_ts_index_caches[ithread][icache]);
        }

        for (unsigned int jthread = 0; jthread < obj->n_threads; ++jthread) {
            free(obj->thread_read_counters[ithread][jthread]);
        }

        free(obj->thread_hare_items_caches[ithread]);
        free(obj->thread_hare_ts_index_caches[ithread]);
    }

    free(obj->thread_tortoise_items);
    free(obj->thread_tortoise_ts_indices);
    free(obj->thread_hare_items_caches);
    free(obj->thread_hare_ts_index_caches);

    for (unsigned int ithread = 0; ithread < obj->n_threads; ++ithread) {
        free(obj->thread_read_counters[ithread]);
        free(obj->thread_write_counters[ithread]);
    }

    free(obj->thread_read_counters);
    free(obj->thread_write_counters);

    for (unsigned int irand = 0; irand < obj->n_rand_items; ++irand) {
        free(obj->random_tG_add_skG[irand]);
        free(obj->random_ts[irand]);
    }

    free(obj->random_tG_add_skG);
    free(obj->random_ts);

    for (unsigned int ithread = 0; ithread < obj->n_threads; ++ithread) {
        free(obj->thread_result_tortoise_ts_indices[ithread]);
        free(obj->thread_result_hare_ts_indices[ithread]);
    }
    
    free(obj->thread_result_tortoise_ts_indices);
    free(obj->thread_result_hare_ts_indices);

    free(obj->founds);

    #ifdef DLOG_VERBOSE
        for (unsigned int ithread = 0; ithread < obj->n_threads; ++ithread) {
            mpz_clear(obj->thread_cache_hit_counters[ithread]);
            mpz_clear(obj->thread_cache_miss_counters[ithread]);
            mpz_clear(obj->thread_cache_possible_misread_counters[ithread]);
        }
        free(obj->thread_cache_hit_counters);
        free(obj->thread_cache_miss_counters);
        free(obj->thread_cache_possible_misread_counters);
    #endif
}

void dlog_print_cache_performance_report(
    dlog_obj obj
)
{
    #ifdef DLOG_VERBOSE
        mpz_t total;
        mpz_init(total);

        printf("\n*********************** CACHE REPORT ***********************\n");
        printf("[debug] CACHE hit all threads: ");
        mpz_set_ui(total, 0);
        for (unsigned int i = 0; i < obj->n_threads; ++i)
            mpz_add(total, total, obj->thread_cache_hit_counters[i]);
        mpz_out_str(stdout, 10, total);
        printf("\n");

        for (unsigned int i = 0; i < obj->n_threads; ++i) {
            printf("[debug]    -> thread %d: ", i);
            mpz_out_str(stdout, 10, obj->thread_cache_hit_counters[i]);
            printf("\n");
        }

        printf("[debug] CACHE miss all threads: ");
        mpz_set_ui(total, 0);
        for (unsigned int i = 0; i < obj->n_threads; ++i)
            mpz_add(total, total, obj->thread_cache_miss_counters[i]);
        mpz_out_str(stdout, 10, total);
        printf("\n");

        for (unsigned int i = 0; i < obj->n_threads; ++i) {
            printf("[debug]    -> thread %d: ", i);
            mpz_out_str(stdout, 10, obj->thread_cache_miss_counters[i]);
            printf("\n");
        }

        printf("[debug] CACHE possible misread all threads: ");
        mpz_set_ui(total, 0);
        for (unsigned int i = 0; i < obj->n_threads; ++i)
            mpz_add(total, total, obj->thread_cache_possible_misread_counters[i]);
        mpz_out_str(stdout, 10, total);
        printf("\n");

        for (unsigned int i = 0; i < obj->n_threads; ++i) {
            printf("[debug]    -> thread %d: ", i);
            mpz_out_str(stdout, 10, obj->thread_cache_possible_misread_counters[i]);
            printf("\n");
        }
        printf("************************************************************\n");
        mpz_clear(total);
    #endif
}

void* __thread__dlog_thread(
    void* vargs
)
{
    // -------------------------------------------------------------------------------------
    //      Setup arguments.
    // -------------------------------------------------------------------------------------
    __args_thread__dlog_thread* args = (__args_thread__dlog_thread*) vargs;
    
    dlog_obj_ptr shared_obj = args->shared_obj;

    unsigned int thread_no     = args->thread_no;
    unsigned int n_threads     = shared_obj->n_threads;
    unsigned int n_cache_items = shared_obj->n_cache_items;
    unsigned int n_rand_items  = shared_obj->n_rand_items;

    mp_size_t item_size_limbs  = shared_obj->item_size_limbs;
    mp_size_t index_size_limbs = shared_obj->index_size_limbs;

    mp_limb_t* field_p  = shared_obj->field_p;
    mp_limb_t* field_P  = shared_obj->field_P;
    mp_limb_t* G_order  = shared_obj->G_order;

    mp_limb_t**  all_tortoise_items        = shared_obj->thread_tortoise_items;
    mp_limb_t**  all_tortoise_ts_indices   = shared_obj->thread_tortoise_ts_indices;
    mp_limb_t*** all_hare_item_caches      = shared_obj->thread_hare_items_caches;
    mp_limb_t*** all_hare_ts_index_caches  = shared_obj->thread_hare_ts_index_caches;

    mp_limb_t*  tortoise_item              = all_tortoise_items[thread_no];
    mp_limb_t*  tortoise_ts_index          = all_tortoise_ts_indices[thread_no];
    mp_limb_t** hare_item_cache            = all_hare_item_caches[thread_no];
    mp_limb_t** hare_ts_index_cache        = all_hare_ts_index_caches[thread_no];

    unsigned long** read_counters      = shared_obj->thread_read_counters[thread_no];
    unsigned long** all_write_counters = shared_obj->thread_write_counters;

    mp_limb_t** random_tG_add_skG = shared_obj->random_tG_add_skG;    
    mp_limb_t** random_ts         = shared_obj->random_ts;

    mp_limb_t* result_tortoise_ts_index = shared_obj->thread_result_tortoise_ts_indices[thread_no];
    mp_limb_t* result_hare_ts_index     = shared_obj->thread_result_hare_ts_indices[thread_no];

    #ifdef DLOG_VERBOSE
        mpz_ptr cache_hit_counter              = shared_obj->thread_cache_hit_counters[thread_no];
        mpz_ptr cache_miss_counter             = shared_obj->thread_cache_miss_counters[thread_no];
        mpz_ptr cache_possible_misread_counter = shared_obj->thread_cache_possible_misread_counters[thread_no];
    #endif

    // -------------------------------------------------------------------------------------
    //      Real calculation.
    //      Do a cycle detection using Brent's algorithm
    //      with Teske's psuedorandom function.
    // -------------------------------------------------------------------------------------
    unsigned long power = 1;
    unsigned long lamda = 1;
    unsigned long icache = 0;

    mp_limb_t* T;
    T = mpn_init_zero(item_size_limbs * 6);

    mp_limb_t* tmp_item;
    mp_limb_t* tmp_ts_index;
    tmp_item = mpn_init_zero(item_size_limbs);
    tmp_ts_index = mpn_init_zero(index_size_limbs * 2);

    while (!shared_obj->overall_found) {
        // ---------------------- updating the tortoise pointer -------------------------
        //                       (every 2^n steps, n increasing)                          
        if (power == lamda) {
            power <<= 1;
            assertf(power != 0, "Oh geez computers are so good nowadays. We've reached the limit of %ld bits, there's nothing we can do now...", sizeof(unsigned long) * 8);
            lamda = 0;

            mpn_copyd(tortoise_item, hare_item_cache[icache], item_size_limbs);
            mpn_copyd(tortoise_ts_index, hare_ts_index_cache[icache], index_size_limbs * 2);
        }

        // ---------------------- updating the hare pointer -------------------------
        unsigned int next_icache = (icache+1) % n_cache_items;
        unsigned int irand = hare_item_cache[icache][0] % n_rand_items;

        // ---------------------- write to a public pointer -------------------------
        //               (this is where the race condition might happen)

        all_write_counters[thread_no][next_icache]++;

        // element <- f(element)
        //     --- and ---
        //      write item
        mpn_montgomery_mulmod_n(
            hare_item_cache[next_icache], 
            hare_item_cache[icache], 
            random_tG_add_skG[irand],

            field_p, field_P,
            item_size_limbs,
            T
        );

        // write updated index.
        mpn_addmod_n(
            &hare_ts_index_cache[next_icache][0],
            &hare_ts_index_cache[icache][0],
            &random_ts[irand][0],
            G_order,
            index_size_limbs
        );

        mpn_addmod_n(
            &hare_ts_index_cache[next_icache][index_size_limbs],
            &hare_ts_index_cache[icache][index_size_limbs],
            &random_ts[irand][index_size_limbs],
            G_order,
            index_size_limbs
        );

        all_write_counters[thread_no][next_icache]++;

        // ------------------- comparing with our hare pointer -------------------
        if (mpn_cmp(tortoise_item, hare_item_cache[next_icache], item_size_limbs) == 0) {
            mpn_copyd(result_tortoise_ts_index, tortoise_ts_index,                index_size_limbs * 2);
            mpn_copyd(result_hare_ts_index,     hare_ts_index_cache[next_icache], index_size_limbs * 2);

            shared_obj->founds[thread_no] = 1;
            shared_obj->overall_found = 1;
            goto dlog_thread_cleanup;
        }

        // ---------- comparing with the other thread's hare pointers -----------
        for (unsigned int ithread = 0; ithread < n_threads; ++ithread) {
            if (ithread == thread_no)
                continue;

            for (unsigned int icache = 0; icache < n_cache_items; ++icache) {
                // If same index as other thread's cache,
                // just move on.
                if (read_counters[ithread][icache] == all_write_counters[ithread][icache]) {
                    #ifdef DLOG_VERBOSE
                        mpz_add_ui(cache_hit_counter, cache_hit_counter, 1);
                    #endif
                    continue;
                }

                /* If the last bit of a write counter
                is 1, that means that the thread is writing,
                we can just skip it. */
                if (all_write_counters[ithread][icache] & 0x1) {
                    continue;
                }

                // Check how many miss happened
                #ifdef DLOG_VERBOSE
                    mpz_add_ui(
                        cache_miss_counter, 
                        cache_miss_counter, 
                        (all_write_counters[ithread][icache] >> 1) - (read_counters[ithread][icache] >> 1) - 1
                    );
                #endif 

                // Set read counter.
                read_counters[ithread][icache] = all_write_counters[ithread][icache];

                // Copy point to avoid the race condition
                // as much as possible.
                mpn_copyd(tmp_ts_index, all_hare_ts_index_caches[ithread][icache], index_size_limbs * 2);
                mpn_copyd(tmp_item, all_hare_item_caches[ithread][icache], item_size_limbs);

                // Check how many misreads might occur.
                #ifdef DLOG_VERBOSE
                    if (read_counters[ithread][icache] != all_write_counters[ithread][icache]) {
                        mpz_add_ui(
                            cache_possible_misread_counter, 
                            cache_possible_misread_counter, 
                            1
                        );
                    }
                #endif 

                // Compare point.
                if (mpn_cmp(tortoise_item, tmp_item, item_size_limbs) == 0) {
                    mpn_copyd(result_tortoise_ts_index, tortoise_ts_index, index_size_limbs * 2);
                    mpn_copyd(result_hare_ts_index,     tmp_ts_index,      index_size_limbs * 2);

                    shared_obj->founds[thread_no] = 1;
                    shared_obj->overall_found = 1;
                    goto dlog_thread_cleanup;
                }
            }
        }

        // ---------------------------- updating lambda --------------------------
        lamda++;
        icache = next_icache;
    }

    // -------------------------------------------------------------------------------------
    //      Cleanup
    // -------------------------------------------------------------------------------------
dlog_thread_cleanup:
    ecc_free_ptemp(T);
    free(tmp_item);
    free(tmp_ts_index);
    free(T);

    return NULL;
}

void dlog_cycle_search(
    dlog_obj obj
)
{
    // -------------------------------------------------------------------------------------
    //      Initialize arguments of __thread__dlog_thread.
    // -------------------------------------------------------------------------------------
    __args_thread__dlog_thread* thread_args = (__args_thread__dlog_thread*) malloc_exit_when_null(sizeof(__args_thread__dlog_thread) * obj->n_threads);
    for (unsigned int ithread = 0; ithread < obj->n_threads; ++ithread) {
        thread_args[ithread].shared_obj = obj;
        thread_args[ithread].thread_no  = ithread;
    }

    // -------------------------------------------------------------------------------------
    //      Generate threads.
    // -------------------------------------------------------------------------------------
    pthread_t* threads = (pthread_t*) malloc_exit_when_null(sizeof(pthread_t) * obj->n_threads);
    int result_code;
    for (unsigned int ithread = 0; ithread < obj->n_threads; ++ithread) {
        result_code = pthread_create(&threads[ithread], NULL, __thread__dlog_thread, (void*)(&thread_args[ithread]));
        if (result_code) {
            printf("[error] oh no! dlog_cycle_search cannot CREATE thread!!!\n");
            exit(-1);
        }
    }

    for (unsigned int ithread = 0; ithread < obj->n_threads; ++ithread) {
        result_code = pthread_join(threads[ithread], NULL);
        if (result_code) {
            printf("[error] oh no! dlog_cycle_search cannot JOIN thread!!!\n");
            exit(-1);
        }
    }

    // -------------------------------------------------------------------------------------
    //      Cleaning up.
    // -------------------------------------------------------------------------------------
    free(thread_args);
    free(threads);
}

void dlog_reset_search(
    dlog_obj obj
)
{
    memset(obj->founds, 0, sizeof(int) * obj->n_threads);
    obj->overall_found = 0;
}

int dlog_get_answer(
    mpz_t p,
    mpz_t k,
    mpz_t G, mpz_t kG,
    mpz_t G_mult_order,

    dlog_obj obj
)
{
    int dlog_status = DLOG_BAD_COLLISION;

    mpz_t t1; mpz_init(t1);
    mpz_t s1; mpz_init(s1);
    mpz_t t2; mpz_init(t2);
    mpz_t s2; mpz_init(s2);

    mpz_t tmp; mpz_init(tmp);

    for (unsigned int ithread = 0; ithread < obj->n_threads; ++ithread) {
        if (obj->founds[ithread]) {
            mpz_set_mpn(t1, &obj->thread_result_tortoise_ts_indices[ithread][0],                     obj->index_size_limbs);
            mpz_set_mpn(s1, &obj->thread_result_tortoise_ts_indices[ithread][obj->index_size_limbs], obj->index_size_limbs);
            mpz_set_mpn(t2, &obj->thread_result_hare_ts_indices[ithread][0],                         obj->index_size_limbs);
            mpz_set_mpn(s2, &obj->thread_result_hare_ts_indices[ithread][obj->index_size_limbs],     obj->index_size_limbs);

            // tortoise = hare
            if (mpz_cmp(t1, t2) != 0) {
                mpz_sub(k, s2, s1);
                mpz_invert(k, k, G_mult_order);
                mpz_sub(tmp, t1, t2);
                mpz_mul(k, k, tmp);
                mpz_mod(k, k, G_mult_order);

                // if k*G == kG -> found!
                mpz_powm(tmp, G, k, p);
                if (mpz_cmp(tmp, kG) == 0) {
                    dlog_status = DLOG_SUCCESS;
                    break;
                }
            }
        }
    }

    mpz_clear(t1);
    mpz_clear(s1);
    mpz_clear(t2);
    mpz_clear(s2);
    mpz_clear(tmp);
    return dlog_status;
}

int dlog(
    mpz_t p,
    mpz_t k, 
    mpz_t G, mpz_t kG, 
    mpz_t G_mult_order, 

    unsigned int n_threads,
    unsigned int n_cache_items,
    unsigned int n_rand_items
)
{
    #ifdef DLOG_VERBOSE
        printf("[debug] p: \n");
        printf("[debug]    ");
        mpz_out_str(stdout, 10, p);
        printf("\n");
        printf("[debug] G: \n");
        printf("[debug]    ");
        mpz_out_str(stdout, 10, G);
        printf("\n");
        printf("[debug] kG: \n");
        printf("[debug]    ");
        mpz_out_str(stdout, 10, kG);
        printf("\n");
        printf("[debug] G_mult_order = ");
        mpz_out_str(stdout, 10, G_mult_order);
        printf("\n");
        printf("[debug] n_threads = %d\n", n_threads);
        printf("[debug] n_cache_items = %d\n", n_cache_items);
        printf("[debug] n_rand_items = %d\n", n_rand_items);
    #endif

    int dlog_status;

    // -------------------------------------------------------------------------------------
    //      Verify if input is good.
    // -------------------------------------------------------------------------------------
    dlog_status = dlog_validate_input(
                    p, 
                    G, kG, 
                    G_mult_order, 

                    n_threads, 
                    n_cache_items, 
                    n_rand_items
                  );
    if (dlog_status != DLOG_MOVE_TO_NEXT_STEP)
        return dlog_status;

    // -------------------------------------------------------------------------------------
    //      Quickly determine if we can/or can't solve
    //      and provide an answer if possible. The reason
    //      for this is that the Pollard-rho algorithm
    //      is probabilistic rather than deterministic. Doing
    //      this will prevent the algorithm from spiraling
    //      into an infinite loop.
    // -------------------------------------------------------------------------------------
    dlog_status = dlog_fast_solve_if_possible(
                    p, 
                    k, 
                    G, kG, 
                    G_mult_order
                  );
    if (dlog_status != DLOG_MOVE_TO_NEXT_STEP)
        return dlog_status;

    // -------------------------------------------------------------------------------------
    //      Initialize memory for fast computations.
    //      (i hope it's fast...)
    // -------------------------------------------------------------------------------------
    dlog_obj obj;
    dlog_init_dlog_obj(
        obj,

        p,
        G, kG,
        G_mult_order,

        n_threads,
        n_cache_items,
        n_rand_items
    );

    dlog_fill_dlog_obj(
        obj,

        p,
        G, kG,
        G_mult_order,

        n_threads,
        n_cache_items,
        n_rand_items
    );

    #ifdef DLOG_VERBOSE
        printf("[debug] index_size_limbs = %ld\n", obj->index_size_limbs);
        printf("[debug] item_size_limbs = %ld\n", obj->item_size_limbs);
    #endif

    // -------------------------------------------------------------------------------------
    //      Doing dlog()...
    // -------------------------------------------------------------------------------------
    while (1) {
        dlog_cycle_search(obj);
        #ifdef DLOG_VERBOSE
            printf("[debug] Collision found!\n");
        #endif

        if (dlog_get_answer(p, k, G, kG, G_mult_order, obj) == DLOG_SUCCESS) {
            #ifdef DLOG_VERBOSE
                printf("[debug] Finished. Found k = ");
                mpz_out_str(stdout, 10, k);
                printf("\n");
            #endif
            break;
        }

        #ifdef DLOG_VERBOSE
            printf("[debug] The result is not correct!? Running cycle detection again...\n");
        #endif
        dlog_reset_search(obj);
    }

    // -------------------------------------------------------------------------------------
    //      Print out report...
    // -------------------------------------------------------------------------------------
    #ifdef DLOG_VERBOSE
        dlog_print_cache_performance_report(obj);
    #endif

    // -------------------------------------------------------------------------------------
    //      Bye! There's no way we reach
    //      this step and the algorithm fails :D
    // -------------------------------------------------------------------------------------
    dlog_free_dlog_obj(obj);
    return DLOG_SUCCESS;
}