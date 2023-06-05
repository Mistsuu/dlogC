#include "dlog.h"
#include "ecc.h"
#include "ecc_proj.h"
#include "ex_mpn.h"
#include "ex_mpz.h"
#include "const.h"
#include "que2.h"
#include "mem.h"

int dlog_validate_input(
    ecc curve,
    eccpt G, eccpt kG,
    mpz_t G_mult_order,

    unsigned int n_threads,
    unsigned int n_caches,
    unsigned int n_randindices
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

    if (n_caches == 0) {
        #ifdef DLOG_VERBOSE
            printf("[debug] You must have >= 1 item in the cache. Exiting...\n");
        #endif
        return DLOG_BAD_CONFIG;
    }

    if (n_randindices < 2) {
        #ifdef DLOG_VERBOSE
            printf("[debug] You must have >= 2 random elements. Exiting...\n");
        #endif
        return DLOG_BAD_CONFIG;
    }

    // -------------------------------------------------------------------------------------
    //      Check G, kG on curve.
    // -------------------------------------------------------------------------------------
    if (!ecc_verify_pt(curve, G)) {
        #ifdef DLOG_VERBOSE
            printf("[debug] G is not on the curve. Exiting...\n");
        #endif
        return DLOG_POINT_NOT_ON_CURVE;
    }

    if (!ecc_verify_pt(curve, kG)) {
        #ifdef DLOG_VERBOSE
            printf("[debug] kG is not on the curve. Exiting...\n");
        #endif
        return DLOG_POINT_NOT_ON_CURVE;
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
    //      Check Gorder * G = O.
    //      Not checking Gorder * kG = O here since
    //      it's not a pre-requisite for this algorithm
    //      to run.
    // -------------------------------------------------------------------------------------
    eccpt O;
    ecc_init_pt(O);
    ecc_mul_noverify(curve, O, G, G_mult_order);
    if (mpz_cmp_ui(O->z, 0) != 0) {
        #ifdef DLOG_VERBOSE
            printf("[debug] G_mult_order * G != infinity. Exiting...\n");
        #endif
        ecc_free_pt(O);
        return DLOG_FAULTY_POINT_ORDER;
    }

    ecc_free_pt(O);
    return DLOG_MOVE_TO_NEXT_STEP;
}

int dlog_fast_solve_if_possible(
    ecc curve,
    mpz_t k, 
    eccpt G, eccpt kG,
    mpz_t G_mult_order
)
{
    // -------------------------------------------------------------------------------------
    //      Case 0: kG = O
    // -------------------------------------------------------------------------------------
    if (mpz_cmp_ui(kG->z, 0) == 0) {
            #ifdef DLOG_VERBOSE
            printf("[debug] kG == O. Automatically set k = 0.\n");
        #endif
        mpz_set_ui(k, 0);
        return DLOG_SUCCESS;
    }

    // -------------------------------------------------------------------------------------
    //      Case 1: G == O and kG != O
    // -------------------------------------------------------------------------------------
    if (mpz_cmp_ui(G->z, 0) == 0) {
        #ifdef DLOG_VERBOSE
            printf("[debug] G == O and kG != O. No solution exists.\n");
        #endif
        return DLOG_NOT_FOUND_DLOG;
    }

    // -------------------------------------------------------------------------------------
    //      Case 2: kG.order() > G.order()
    // -------------------------------------------------------------------------------------
    eccpt O;
    ecc_init_pt(O);
    ecc_mul_noverify(curve, O, kG, G_mult_order);
    if (mpz_cmp_ui(O->z, 0) != 0) {
        #ifdef DLOG_VERBOSE
            printf("[debug] n*G == infinity but n*kG != infinity. No solution exists.\n");
        #endif
        ecc_free_pt(O);
        return DLOG_NOT_FOUND_DLOG;
    }

    // -------------------------------------------------------------------------------------
    //      Case 3: Using Weil Pairing to 
    //              check if there's solution?
    // -------------------------------------------------------------------------------------
    mpz_t E;
    mpz_init(E);
    ecc_weil_pairing_noverify(curve, E, G, kG, G_mult_order);
    if (mpz_cmp_ui(E, 1) != 0) {
        #ifdef DLOG_VERBOSE
            printf("[debug] weil_pairing(G, kG) != 1. kG is not a multiple of G.\n");
        #endif
        mpz_clear(E);
        ecc_free_pt(O);
        return DLOG_NOT_FOUND_DLOG;
    }

    mpz_clear(E);
    ecc_free_pt(O);
    return DLOG_MOVE_TO_NEXT_STEP;
}

void dlog_init_dlog_obj(
    dlog_obj obj,

    ecc curve,
    eccpt G, eccpt kG,
    mpz_t G_mult_order,

    unsigned int n_threads,
    unsigned int n_caches,
    unsigned int n_randindices
)
{
    obj->n_threads = n_threads;
    obj->n_caches = n_caches;
    obj->n_randindices = n_randindices;
    obj->item_size_limbs  = mpz_size(curve->p);
    obj->index_size_limbs = mpz_size(G_mult_order);

    obj->curve_aR = mpn_init_zero(obj->item_size_limbs);
    obj->curve_bR = mpn_init_zero(obj->item_size_limbs);
    obj->curve_p  = mpn_init_zero(obj->item_size_limbs);
    obj->curve_P  = mpn_init_zero(obj->item_size_limbs);
    obj->G_order  = mpn_init_zero(obj->index_size_limbs);

    obj->thread_tortoise_item_caches = (mp_limb_t***)malloc_exit_when_null(sizeof(mp_limb_t**) * n_threads);
    obj->thread_hare_item_caches = (mp_limb_t***)malloc_exit_when_null(sizeof(mp_limb_t**) * n_threads);
    obj->thread_tortoise_index_caches = (mp_limb_t***)malloc_exit_when_null(sizeof(mp_limb_t**) * n_threads);
    obj->thread_hare_index_caches = (mp_limb_t***)malloc_exit_when_null(sizeof(mp_limb_t**) * n_threads);
    obj->thread_read_counters = (unsigned long***)malloc_exit_when_null(sizeof(unsigned long**) * n_threads);
    obj->thread_write_counters = (unsigned long**)malloc_exit_when_null(sizeof(unsigned long*) * n_threads);

    for (unsigned int i = 0; i < n_threads; ++i) {
        obj->thread_tortoise_item_caches[i] = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_caches);
        obj->thread_hare_item_caches[i] = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_caches);
        obj->thread_tortoise_index_caches[i] = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_caches);
        obj->thread_hare_index_caches[i] = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_caches);
        obj->thread_read_counters[i] = (unsigned long**)malloc_exit_when_null(sizeof(unsigned long*) * n_threads);
        obj->thread_write_counters[i] = (unsigned long*)malloc_exit_when_null(sizeof(unsigned long) * n_caches);

        for (unsigned int j = 0; j < n_caches; ++j) {
            obj->thread_tortoise_item_caches[i][j] = mpn_init_zero(obj->item_size_limbs * 3);
            obj->thread_hare_item_caches[i][j] = mpn_init_zero(obj->item_size_limbs * 3);
            obj->thread_tortoise_index_caches[i][j] = mpn_init_zero(obj->index_size_limbs * 2);
            obj->thread_hare_index_caches[i][j] = mpn_init_zero(obj->index_size_limbs * 2);
        }

        for (unsigned int j = 0; j < n_threads; ++j) {
            obj->thread_read_counters[i][j] = (unsigned long*)malloc_exit_when_null(sizeof(unsigned long) * n_caches);
        }
    }

    obj->random_tG_add_skG = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_randindices);
    obj->random_ts         = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_randindices);
    for (unsigned int i = 0; i < n_randindices; ++i) {
        obj->random_tG_add_skG[i] = mpn_init_zero(obj->item_size_limbs * 2);
        obj->random_ts[i] = mpn_init_zero(obj->index_size_limbs * 2);
    }

    obj->thread_result_tortoise_items   = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_threads);
    obj->thread_result_hare_items       = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_threads);
    obj->thread_result_tortoise_indices = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_threads);
    obj->thread_result_hare_indices     = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_threads);
    for (unsigned int i = 0; i < n_threads; ++i) {
        obj->thread_result_tortoise_items[i]   = mpn_init_zero(obj->item_size_limbs * 3);
        obj->thread_result_hare_items[i]       = mpn_init_zero(obj->item_size_limbs * 3);
        obj->thread_result_tortoise_indices[i] = mpn_init_zero(obj->index_size_limbs * 2);
        obj->thread_result_hare_indices[i]     = mpn_init_zero(obj->index_size_limbs * 2);
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

    ecc curve,
    eccpt G, eccpt kG,
    mpz_t G_mult_order,

    unsigned int n_threads,
    unsigned int n_caches,
    unsigned int n_randindices
)
{
    mpz_t mpz_R;
    mpz_init_set_ui(mpz_R, 1);
    mpz_mul_2exp(mpz_R, mpz_R, obj->item_size_limbs * mp_bits_per_limb);

    // -------------------------------------------------------------------------------------
    //      We have to convert curve's a and curve's b
    //      to Montgomery form.
    // -------------------------------------------------------------------------------------
    mpz_t mpz_aR;
    mpz_t mpz_bR;
    mpz_t mpz_curve_P;
    mpz_init(mpz_aR);
    mpz_init(mpz_bR);
    mpz_init(mpz_curve_P);
    mpz_mul(mpz_aR, curve->a, mpz_R);
    mpz_mul(mpz_bR, curve->b, mpz_R);
    mpz_mod(mpz_aR, mpz_aR, curve->p);
    mpz_mod(mpz_bR, mpz_bR, curve->p);
    mpz_invert(mpz_curve_P, curve->p, mpz_R);

    mpn_cpyz(obj->curve_aR, mpz_aR,       obj->item_size_limbs);
    mpn_cpyz(obj->curve_bR, mpz_bR,       obj->item_size_limbs);
    mpn_cpyz(obj->curve_p,  curve->p,     obj->item_size_limbs);
    mpn_cpyz(obj->curve_P,  mpz_curve_P,  obj->item_size_limbs);
    mpn_cpyz(obj->G_order,  G_mult_order, obj->index_size_limbs);

    // -------------------------------------------------------------------------------------
    //      Initialize fixed random points
    // -------------------------------------------------------------------------------------
    eccpt tG;
    eccpt skG;
    eccpt tG_add_skG;
    mpz_t t;
    mpz_t s;
    mpz_init(t);
    mpz_init(s);
    ecc_init_pt(tG);
    ecc_init_pt(skG);
    ecc_init_pt(tG_add_skG);
    
    for (unsigned int i = 0; i < n_randindices; ++i) {
        mpz_dev_urandomm(t, G_mult_order);
        mpz_dev_urandomm(s, G_mult_order);

        ecc_mul_noverify(curve, tG, G, t);
        ecc_mul_noverify(curve, skG, kG, s);
        ecc_add_noverify(curve, tG_add_skG, tG, skG);

        // Convert coordinate to Montgomery form.
        mpz_mul(tG_add_skG->x, tG_add_skG->x, mpz_R);
        mpz_mod(tG_add_skG->x, tG_add_skG->x, curve->p);
        mpz_mul(tG_add_skG->y, tG_add_skG->y, mpz_R);
        mpz_mod(tG_add_skG->y, tG_add_skG->y, curve->p);

        mpn_cpyz( obj->random_ts[i],                        t, obj->index_size_limbs);
        mpn_cpyz(&obj->random_ts[i][obj->index_size_limbs], s, obj->index_size_limbs);
        mpn_cpyz( obj->random_tG_add_skG[i],                       tG_add_skG->x, obj->item_size_limbs);
        mpn_cpyz(&obj->random_tG_add_skG[i][obj->item_size_limbs], tG_add_skG->y, obj->item_size_limbs);
    }

    // -------------------------------------------------------------------------------------
    //      Initialize cache values
    // -------------------------------------------------------------------------------------
    mpz_t mpz_1;
    mpz_init_set_ui(mpz_1, 1);

    for (unsigned int i = 0; i < n_threads; ++i) {
        for (unsigned int j = 0; j < n_caches; ++j) {
            mpz_dev_urandomm(t, G_mult_order);
            mpz_dev_urandomm(s, G_mult_order);

            ecc_mul_noverify(curve, tG, G, t);
            ecc_mul_noverify(curve, skG, kG, s);
            ecc_add_noverify(curve, tG_add_skG, tG, skG);

            mpn_cpyz( obj->thread_tortoise_item_caches[i][j],                         tG_add_skG->x, obj->item_size_limbs);
            mpn_cpyz(&obj->thread_tortoise_item_caches[i][j][obj->item_size_limbs],   tG_add_skG->y, obj->item_size_limbs);
            mpn_cpyz(&obj->thread_tortoise_item_caches[i][j][obj->item_size_limbs*2], mpz_1,         obj->item_size_limbs);

            mpn_cpyz( obj->thread_hare_item_caches[i][j],                         tG_add_skG->x, obj->item_size_limbs);
            mpn_cpyz(&obj->thread_hare_item_caches[i][j][obj->item_size_limbs],   tG_add_skG->y, obj->item_size_limbs);
            mpn_cpyz(&obj->thread_hare_item_caches[i][j][obj->item_size_limbs*2], mpz_1,         obj->item_size_limbs);

            mpn_cpyz( obj->thread_tortoise_index_caches[i][j],                        t, obj->index_size_limbs);
            mpn_cpyz(&obj->thread_tortoise_index_caches[i][j][obj->index_size_limbs], s, obj->index_size_limbs);

            mpn_cpyz( obj->thread_hare_index_caches[i][j],                        t, obj->index_size_limbs);
            mpn_cpyz(&obj->thread_hare_index_caches[i][j][obj->index_size_limbs], s, obj->index_size_limbs);
        }
    }

    // -------------------------------------------------------------------------------------
    //      Initialize read/write counters
    // -------------------------------------------------------------------------------------
    for (unsigned int i = 0; i < n_threads; ++i) {
        for (unsigned int j = 0; j < n_caches; ++j) {
            obj->thread_write_counters[i][j] = 0;
        }

        for (unsigned int j = 0; j < n_threads; ++j) {
            for (unsigned int k = 0; k < n_caches; ++k) {
                obj->thread_read_counters[i][j][k] = 0;
            }
        }
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
        for (unsigned int i = 0; i < n_threads; ++i) {
            mpz_set_ui(obj->thread_cache_hit_counters[i], 0);
            mpz_set_ui(obj->thread_cache_miss_counters[i], 0);
            mpz_set_ui(obj->thread_cache_possible_misread_counters[i], 0);
        }
    #endif

    // -------------------------------------------------------------------------------------
    //      Free stuffs
    // -------------------------------------------------------------------------------------
    mpz_clear(mpz_R);
    mpz_clear(mpz_aR);
    mpz_clear(mpz_bR);
    mpz_clear(mpz_curve_P);
    mpz_clear(mpz_1);
    mpz_clear(t);
    mpz_clear(s);
    ecc_free_pt(tG);
    ecc_free_pt(skG);
    ecc_free_pt(tG_add_skG);
}

void dlog_free_dlog_obj(
    dlog_obj obj
)
{
    free(obj->curve_aR);
    free(obj->curve_bR);
    free(obj->curve_p);
    free(obj->curve_P);
    free(obj->G_order);

    for (unsigned int i = 0; i < obj->n_threads; ++i) {
        for (unsigned int j = 0; j < obj->n_caches; ++j) {
            free(obj->thread_tortoise_item_caches[i][j]);
            free(obj->thread_hare_item_caches[i][j]);
            free(obj->thread_tortoise_index_caches[i][j]);
            free(obj->thread_hare_index_caches[i][j]);
        }

        for (unsigned int j = 0; j < obj->n_threads; ++j) {
            free(obj->thread_read_counters[i][j]);
        }

        free(obj->thread_tortoise_item_caches[i]);
        free(obj->thread_hare_item_caches[i]);
        free(obj->thread_tortoise_index_caches[i]);
        free(obj->thread_hare_index_caches[i]);
        free(obj->thread_read_counters[i]);
        free(obj->thread_write_counters[i]);
    }

    free(obj->thread_tortoise_item_caches);
    free(obj->thread_hare_item_caches);
    free(obj->thread_tortoise_index_caches);
    free(obj->thread_hare_index_caches);
    free(obj->thread_read_counters);
    free(obj->thread_write_counters);

    for (unsigned int i = 0; i < obj->n_randindices; ++i) {
        free(obj->random_tG_add_skG[i]);
        free(obj->random_ts[i]);
    }

    free(obj->random_tG_add_skG);
    free(obj->random_ts);

    for (unsigned int i = 0; i < obj->n_threads; ++i) {
        free(obj->thread_result_tortoise_items[i]);
        free(obj->thread_result_hare_items[i]);
        free(obj->thread_result_tortoise_indices[i]);
        free(obj->thread_result_hare_indices[i]);
    }
    
    free(obj->thread_result_tortoise_items);
    free(obj->thread_result_hare_items);
    free(obj->thread_result_tortoise_indices);
    free(obj->thread_result_hare_indices);

    free(obj->founds);

    #ifdef DLOG_VERBOSE
        for (unsigned int i = 0; i < obj->n_threads; ++i) {
            mpz_clear(obj->thread_cache_hit_counters[i]);
            mpz_clear(obj->thread_cache_miss_counters[i]);
            mpz_clear(obj->thread_cache_possible_misread_counters[i]);
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

}

int dlog_main_process(
    dlog_obj obj
)
{
    // -------------------------------------------------------------------------------------
    //      Initialize arguments of __thread__dlog_thread.
    // -------------------------------------------------------------------------------------
    __args_thread__dlog_thread* thread_args = (__args_thread__dlog_thread*) malloc_exit_when_null(sizeof(__args_thread__dlog_thread) * obj->n_threads);
    for (unsigned int i = 0; i < obj->n_threads; ++i) {
        thread_args[i].shared_obj = obj;
        thread_args[i].thread_no  = i;
    }

    // -------------------------------------------------------------------------------------
    //      Generate threads.
    // -------------------------------------------------------------------------------------
    pthread_t* threads = (pthread_t*) malloc_exit_when_null(sizeof(pthread_t) * obj->n_threads);
    int result_code;
    for (unsigned int i = 0; i < obj->n_threads; ++i) {
        result_code = pthread_create(&threads[i], NULL, __thread__dlog_thread, (void*)(&thread_args[i]));
        if (result_code) {
            printf("[error] oh no! dlog_main_process cannot CREATE thread!!!\n");
            exit(-1);
        }
    }

    for (unsigned int i = 0; i < obj->n_threads; ++i) {
        result_code = pthread_join(threads[i], NULL);
        if (result_code) {
            printf("[error] oh no! dlog_main_process cannot JOIN thread!!!\n");
            exit(-1);
        }
    }

    // -------------------------------------------------------------------------------------
    //      Cleaning up.
    // -------------------------------------------------------------------------------------
    free(thread_args);
    free(threads);

    // todo: return a variable `dlog_status`.
    return DLOG_SUCCESS;
}

int dlog2(
    ecc curve, 
    mpz_t k, 
    eccpt G, eccpt kG, 
    mpz_t G_mult_order, 

    unsigned int n_threads,
    unsigned int n_caches,
    unsigned int n_randindices
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
        printf("[debug] G_mult_order = ");
        mpz_out_str(stdout, 10, G_mult_order);
        printf("\n");
        printf("[debug] n_threads = %d\n", n_threads);
        printf("[debug] n_caches = %d\n", n_caches);
        printf("[debug] n_randindices = %d\n", n_randindices);
    #endif

    int dlog_status;

    // -------------------------------------------------------------------------------------
    //      Verify if input is good.
    // -------------------------------------------------------------------------------------
    dlog_status = dlog_validate_input(
                    curve, 
                    G, kG, 
                    G_mult_order, 

                    n_threads, 
                    n_caches, 
                    n_randindices
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
                    curve, 
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

        curve,
        G, kG,
        G_mult_order,

        n_threads,
        n_caches,
        n_randindices
    );

    dlog_fill_dlog_obj(
        obj,

        curve,
        G, kG,
        G_mult_order,

        n_threads,
        n_caches,
        n_randindices
    );

    #ifdef DLOG_VERBOSE
        printf("[debug] index_size_limbs = %ld\n", obj->index_size_limbs);
        printf("[debug] item_size_limbs = %ld\n", obj->item_size_limbs);
    #endif

    // -------------------------------------------------------------------------------------
    //      Doing dlog()...
    // -------------------------------------------------------------------------------------
    while (dlog_main_process(obj) == DLOG_BAD_COLLISION) {
        #ifdef DLOG_VERBOSE
            printf("[debug] It seems like we have found collision but the result is not correct! Doing dlog_main_process() again...\n");
        #endif
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