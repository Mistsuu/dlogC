#ifndef DLOG_H
#define DLOG_H

#include <gmp.h>
#include <ex_assert.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>

#ifdef DLOG_VERBOSE
    #include <sys/time.h>
#endif


typedef struct dlog_obj_struct
{   
    unsigned int n_threads;
    unsigned int n_cache_items;
    unsigned int n_rand_items;
    mp_size_t item_size_limbs;
    mp_size_t index_size_limbs;

    /* field */
    mp_limb_t* field_p;
    mp_limb_t* field_P;
    mp_limb_t* G_order;

    /* sorry for all the 3 stars pointers... */
    mp_limb_t**       thread_tortoise_items;            // n threads, each thread has 1 value.
    mp_limb_t***      thread_hare_items_caches;         // n threads, each thread has m cache values.
    mp_limb_t**       thread_tortoise_ts_indices;       // n threads, each thread has 1 value.
    mp_limb_t***      thread_hare_ts_index_caches;      // n threads, each thread has m cache indices.

    unsigned long*    thread_write_index;

    mp_limb_t** random_tG_add_skG;                      // r random points of t*G + s*kG (using X, Y coordinate in Montgomery Form.)
    mp_limb_t** random_ts;                              // r random multipliers t, s

    /* result storage */
    mp_limb_t** thread_result_tortoise_ts_indices;
    mp_limb_t** thread_result_hare_ts_indices;

    /* result indicators */
    int* founds;
    int overall_found;

} __dlog_obj_struct;

typedef __dlog_obj_struct  dlog_obj[1];
typedef __dlog_obj_struct* dlog_obj_ptr;

typedef struct {
    dlog_obj_ptr shared_obj;
    unsigned int thread_no;
} __args_thread__dlog_thread;

int dlog_validate_input(
    mpz_t p,
    mpz_t G, mpz_t kG,
    mpz_t G_mult_order,

    unsigned int n_threads,
    unsigned int n_cache_items,
    unsigned int n_rand_items
);

int dlog_fast_solve_if_possible(
    mpz_t p,
    mpz_t k, 
    mpz_t G, mpz_t kG,
    mpz_t G_mult_order
);

void dlog_init_dlog_obj(
    dlog_obj obj,

    mpz_t p,
    mpz_t G, mpz_t kG,
    mpz_t G_mult_order,

    unsigned int n_threads,
    unsigned int n_cache_items,
    unsigned int n_rand_items
);

void dlog_fill_dlog_obj(
    dlog_obj obj,

    mpz_t p,
    mpz_t G, mpz_t kG,
    mpz_t G_mult_order,

    unsigned int n_threads,
    unsigned int n_cache_items,
    unsigned int n_rand_items
);

void dlog_free_dlog_obj(
    dlog_obj obj
);

void dlog_cycle_search(
    dlog_obj obj
);

void dlog_reset_search(
    dlog_obj obj
);

int dlog_get_answer(
    mpz_t p,
    mpz_t k,
    mpz_t G, mpz_t kG,
    mpz_t G_mult_order,

    dlog_obj obj
);

int dlog(
    mpz_t p,
    mpz_t k, 
    mpz_t G, mpz_t kG, 
    mpz_t G_mult_order, 

    unsigned int n_threads,
    unsigned int n_cache_items,
    unsigned int n_rand_items
);

#endif