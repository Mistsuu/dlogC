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

    obj->thread_tortoise_X_items = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_threads);
    obj->thread_tortoise_ts_indices = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_threads);
    obj->thread_hare_XYZ_items = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_threads);
    obj->thread_hare_X_items_caches = (mp_limb_t***)malloc_exit_when_null(sizeof(mp_limb_t**) * n_threads);
    obj->thread_hare_ts_index_caches = (mp_limb_t***)malloc_exit_when_null(sizeof(mp_limb_t**) * n_threads);

    for (unsigned int ithread = 0; ithread < n_threads; ++ithread) {
        obj->thread_tortoise_X_items[ithread] = mpn_init_zero(obj->item_size_limbs);
        obj->thread_tortoise_ts_indices[ithread] = mpn_init_zero(obj->index_size_limbs * 2);
        
        obj->thread_hare_XYZ_items[ithread] = mpn_init_zero(obj->item_size_limbs * 3);
        obj->thread_hare_X_items_caches[ithread] = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_caches);
        obj->thread_hare_ts_index_caches[ithread] = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_caches);
        for (unsigned int jthread = 0; jthread < n_caches; ++jthread) {
            obj->thread_hare_X_items_caches[ithread][jthread] = mpn_init_zero(obj->item_size_limbs);
            obj->thread_hare_ts_index_caches[ithread][jthread] = mpn_init_zero(obj->index_size_limbs * 2);
        }
    }

    obj->thread_read_counters = (unsigned long***)malloc_exit_when_null(sizeof(unsigned long**) * n_threads);
    obj->thread_write_counters = (unsigned long**)malloc_exit_when_null(sizeof(unsigned long*) * n_threads);

    for (unsigned int ithread = 0; ithread < n_threads; ++ithread) {
        obj->thread_read_counters[ithread] = (unsigned long**)malloc_exit_when_null(sizeof(unsigned long*) * n_threads);
        obj->thread_write_counters[ithread] = (unsigned long*)malloc_exit_when_null(sizeof(unsigned long) * n_caches);
        for (unsigned int jthread = 0; jthread < n_threads; ++jthread) {
            obj->thread_read_counters[ithread][jthread] = (unsigned long*)malloc_exit_when_null(sizeof(unsigned long) * n_caches);
        }
    }

    obj->random_tG_add_skG = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_randindices);
    obj->random_ts         = (mp_limb_t**)malloc_exit_when_null(sizeof(mp_limb_t*) * n_randindices);
    for (unsigned int irand = 0; irand < n_randindices; ++irand) {
        obj->random_tG_add_skG[irand] = mpn_init_zero(obj->item_size_limbs * 2);
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
    mpz_sub(mpz_curve_P, mpz_R, mpz_curve_P);
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
    
    for (unsigned int irand = 0; irand < n_randindices; ++irand) {
        do {
            mpz_dev_urandomm(t, G_mult_order);
            mpz_dev_urandomm(s, G_mult_order);

            ecc_mul_noverify(curve, tG, G, t);
            ecc_mul_noverify(curve, skG, kG, s);
            ecc_add_noverify(curve, tG_add_skG, tG, skG);
        } while (mpz_cmp_ui(tG_add_skG->z, 0) == 0);

        // Convert coordinate to Montgomery form.
        mpz_mul(tG_add_skG->x, tG_add_skG->x, mpz_R);
        mpz_mod(tG_add_skG->x, tG_add_skG->x, curve->p);
        mpz_mul(tG_add_skG->y, tG_add_skG->y, mpz_R);
        mpz_mod(tG_add_skG->y, tG_add_skG->y, curve->p);

        mpn_cpyz( obj->random_ts[irand],                        t, obj->index_size_limbs);
        mpn_cpyz(&obj->random_ts[irand][obj->index_size_limbs], s, obj->index_size_limbs);
        mpn_cpyz( obj->random_tG_add_skG[irand],                       tG_add_skG->x, obj->item_size_limbs);
        mpn_cpyz(&obj->random_tG_add_skG[irand][obj->item_size_limbs], tG_add_skG->y, obj->item_size_limbs);
    }

    // -------------------------------------------------------------------------------------
    //      Initialize cache values
    // -------------------------------------------------------------------------------------
    mpz_t mpz_1;
    mpz_init_set_ui(mpz_1, 1);

    for (unsigned int ithread = 0; ithread < n_threads; ++ithread) {
        do {
            mpz_dev_urandomm(t, G_mult_order);
            mpz_dev_urandomm(s, G_mult_order);

            ecc_mul_noverify(curve, tG, G, t);
            ecc_mul_noverify(curve, skG, kG, s);
            ecc_add_noverify(curve, tG_add_skG, tG, skG);
        } while (mpz_cmp_ui(tG_add_skG->z, 0) == 0);

        mpn_cpyz( obj->thread_tortoise_X_items[ithread], tG_add_skG->x, obj->item_size_limbs);
        mpn_cpyz( obj->thread_tortoise_ts_indices[ithread],                        t, obj->index_size_limbs);
        mpn_cpyz(&obj->thread_tortoise_ts_indices[ithread][obj->index_size_limbs], s, obj->index_size_limbs);

        mpn_cpyz( obj->thread_hare_XYZ_items[ithread],                         tG_add_skG->x, obj->item_size_limbs);
        mpn_cpyz(&obj->thread_hare_XYZ_items[ithread][obj->item_size_limbs],   tG_add_skG->y, obj->item_size_limbs);
        mpn_cpyz(&obj->thread_hare_XYZ_items[ithread][obj->item_size_limbs*2], mpz_1,         obj->item_size_limbs);
        for (unsigned int icache = 0; icache < n_caches; ++icache) {
            mpn_cpyz( obj->thread_hare_X_items_caches[ithread][icache], tG_add_skG->x, obj->item_size_limbs);
            mpn_cpyz( obj->thread_hare_ts_index_caches[ithread][icache],                        t, obj->index_size_limbs);
            mpn_cpyz(&obj->thread_hare_ts_index_caches[ithread][icache][obj->index_size_limbs], s, obj->index_size_limbs);
        }
    }

    // -------------------------------------------------------------------------------------
    //      Initialize read/write counters
    // -------------------------------------------------------------------------------------
    for (unsigned int ithread = 0; ithread < n_threads; ++ithread) {
        memset(obj->thread_write_counters[ithread], 0, n_caches * sizeof(unsigned long));
        for (unsigned int jthread = 0; jthread < n_threads; ++jthread)
            memset(obj->thread_read_counters[ithread][jthread], 0, n_caches * sizeof(unsigned long));
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

    for (unsigned int ithread = 0; ithread < obj->n_threads; ++ithread) {
        free(obj->thread_tortoise_X_items[ithread]);
        free(obj->thread_tortoise_ts_indices[ithread]);

        free(obj->thread_hare_XYZ_items[ithread]);
        for (unsigned int icache = 0; icache < obj->n_caches; ++icache) {
            free(obj->thread_hare_X_items_caches[ithread][icache]);
            free(obj->thread_hare_ts_index_caches[ithread][icache]);
        }

        for (unsigned int jthread = 0; jthread < obj->n_threads; ++jthread) {
            free(obj->thread_read_counters[ithread][jthread]);
        }

        free(obj->thread_hare_X_items_caches[ithread]);
        free(obj->thread_hare_ts_index_caches[ithread]);
    }

    free(obj->thread_tortoise_X_items);
    free(obj->thread_tortoise_ts_indices);
    free(obj->thread_hare_XYZ_items);
    free(obj->thread_hare_X_items_caches);
    free(obj->thread_hare_ts_index_caches);

    for (unsigned int ithread = 0; ithread < obj->n_threads; ++ithread) {
        free(obj->thread_read_counters[ithread]);
        free(obj->thread_write_counters[ithread]);
    }

    free(obj->thread_read_counters);
    free(obj->thread_write_counters);

    for (unsigned int irand = 0; irand < obj->n_randindices; ++irand) {
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
    unsigned int n_caches      = shared_obj->n_caches;
    unsigned int n_randindices = shared_obj->n_randindices;

    mp_size_t item_size_limbs  = shared_obj->item_size_limbs;
    mp_size_t index_size_limbs = shared_obj->index_size_limbs;

    mp_limb_t* curve_aR = shared_obj->curve_aR;
    mp_limb_t* curve_bR = shared_obj->curve_bR;
    mp_limb_t* curve_p  = shared_obj->curve_p;
    mp_limb_t* curve_P  = shared_obj->curve_P;
    mp_limb_t* G_order  = shared_obj->G_order;

    mp_limb_t**  all_tortoise_X_items      = shared_obj->thread_tortoise_X_items;
    mp_limb_t**  all_tortoise_ts_indices   = shared_obj->thread_tortoise_ts_indices;
    mp_limb_t**  all_hare_XYZ_items        = shared_obj->thread_hare_XYZ_items;
    mp_limb_t*** all_hare_X_item_caches    = shared_obj->thread_hare_X_items_caches;
    mp_limb_t*** all_hare_ts_index_caches  = shared_obj->thread_hare_ts_index_caches;

    mp_limb_t*  tortoise_X_item            = all_tortoise_X_items[thread_no];
    mp_limb_t*  tortoise_ts_index          = all_tortoise_ts_indices[thread_no];
    mp_limb_t*  hare_XYZ_item              = all_hare_XYZ_items[thread_no];
    mp_limb_t** hare_X_item_cache          = all_hare_X_item_caches[thread_no];
    mp_limb_t** hare_ts_index_cache        = all_hare_ts_index_caches[thread_no];

    unsigned long** read_counters               = shared_obj->thread_read_counters[thread_no];
    unsigned long** all_write_counters          = shared_obj->thread_write_counters;

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

    ecc_ptemp T;
    ecc_init_ptemp(T, item_size_limbs);

    mp_limb_t* tmp_X_item;
    mp_limb_t* tmp_ts_index;
    tmp_X_item = mpn_init_zero(item_size_limbs * 3);
    tmp_ts_index = mpn_init_zero(index_size_limbs * 2);

    while (!shared_obj->overall_found) {
        // ---------------------- updating the tortoise pointer -------------------------
        //                       (every 2^n steps, n increasing)                          
        if (power == lamda) {
            power <<= 1;
            assertf(power != 0, "Oh geez computers are so good nowadays. We've reached the limit of %ld bits, there's nothing we can do now...", sizeof(unsigned long) * 8);
            lamda = 0;

            mpn_copyd(tortoise_X_item,   hare_X_item_cache[icache],   item_size_limbs);
            mpn_copyd(tortoise_ts_index, hare_ts_index_cache[icache], index_size_limbs * 2);
        }

        // ---------------------- updating the hare pointer -------------------------
        unsigned int next_icache = (icache+1) % n_caches;
        unsigned int irand = hare_X_item_cache[icache][0] % n_randindices; // This should somehow uniquely map x/z -> x.

        // element <- f(element)
        if (!ecc_peq(
                &hare_XYZ_item[0],
                &hare_XYZ_item[item_size_limbs],
                &hare_XYZ_item[item_size_limbs*2],

                &random_tG_add_skG[irand][0],
                &random_tG_add_skG[irand][item_size_limbs],

                curve_p,
                curve_P,
                item_size_limbs,

                T
            ))
                ecc_padd(
                    &hare_XYZ_item[0],
                    &hare_XYZ_item[item_size_limbs],
                    &hare_XYZ_item[item_size_limbs*2],

                    &hare_XYZ_item[0],
                    &hare_XYZ_item[item_size_limbs],
                    &hare_XYZ_item[item_size_limbs*2],

                    &random_tG_add_skG[irand][0],
                    &random_tG_add_skG[irand][item_size_limbs],
                    
                    curve_p,
                    curve_P,
                    item_size_limbs,

                    T
                );
        else
                ecc_pdbl(
                    &hare_XYZ_item[0],
                    &hare_XYZ_item[item_size_limbs],
                    &hare_XYZ_item[item_size_limbs*2],

                    &hare_XYZ_item[0],
                    &hare_XYZ_item[item_size_limbs],
                    &hare_XYZ_item[item_size_limbs*2],

                    curve_aR,
                    curve_p,
                    curve_P,
                    item_size_limbs,

                    T
                );


        // ---------------------- write to a public pointer -------------------------
        //               (this is where the race condition might happen)

        all_write_counters[thread_no][next_icache]++;

        // write normalized X/Z coordinate
        ecc_pxz_to_X(
            hare_X_item_cache[next_icache], 
            
            &hare_XYZ_item[0],
            &hare_XYZ_item[item_size_limbs*2],

            curve_p,
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
        if (mpn_cmp(tortoise_X_item, hare_X_item_cache[next_icache], item_size_limbs) == 0) {
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

            for (unsigned int icache = 0; icache < n_caches; ++icache) {
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
                mpn_copyd(tmp_X_item, all_hare_X_item_caches[ithread][icache], item_size_limbs);

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
                if (mpn_cmp(tortoise_X_item, tmp_X_item, item_size_limbs) == 0) {
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
    free(tmp_X_item);
    free(tmp_ts_index);

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
    ecc curve,
    mpz_t k,
    eccpt G, eccpt kG,
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
    eccpt TMP; ecc_init_pt(TMP);

    for (unsigned int ithread = 0; ithread < obj->n_threads; ++ithread) {
        if (obj->founds[ithread]) {
            mpz_set_mpn(t1, &obj->thread_result_tortoise_ts_indices[ithread][0],                     obj->index_size_limbs);
            mpz_set_mpn(s1, &obj->thread_result_tortoise_ts_indices[ithread][obj->index_size_limbs], obj->index_size_limbs);
            mpz_set_mpn(t2, &obj->thread_result_hare_ts_indices[ithread][0],                         obj->index_size_limbs);
            mpz_set_mpn(s2, &obj->thread_result_hare_ts_indices[ithread][obj->index_size_limbs],     obj->index_size_limbs);

            // because we have x(tortoise) = x(hare)
            // we have 2 different routes:
            // tortoise = hare
            if (mpz_cmp(t1, t2) != 0) {
                mpz_sub(k, s2, s1);
                mpz_invert(k, k, G_mult_order);
                mpz_sub(tmp, t1, t2);
                mpz_mul(k, k, tmp);
                mpz_mod(k, k, G_mult_order);

                // if k*G == kG -> found!
                ecc_mul(curve, TMP, G, k);
                if (ecc_eq(curve, TMP, kG) == 1) {
                    dlog_status = DLOG_SUCCESS;
                    break;
                }
            }

            // tortoise = -hare
            mpz_add(tmp, s1, s2);
            mpz_mod(tmp, tmp, G_mult_order);
            if (mpz_cmp_ui(tmp, 0) != 0) {
                mpz_add(k, s2, s1);
                mpz_invert(k, k, G_mult_order);
                mpz_sub(k, G_mult_order, k);
                mpz_add(tmp, t1, t2);
                mpz_mul(k, k, tmp);
                mpz_mod(k, k, G_mult_order);

                // if k*G == kG -> found!
                ecc_mul(curve, TMP, G, k);
                if (ecc_eq(curve, TMP, kG) == 1) {
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

    ecc_free_pt(TMP);
    return dlog_status;
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
    while (1) {
        dlog_cycle_search(obj);
        #ifdef DLOG_VERBOSE
            printf("[debug] Collision found!\n");
        #endif

        if (dlog_get_answer(curve, k, G, kG, G_mult_order, obj) == DLOG_SUCCESS) {
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