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

#include "ecc.h"
#include "ecc_proj.h"

typedef struct dlog_obj_struct
{   
    unsigned int n_threads;
    unsigned int n_caches;
    unsigned int n_randindices;
    mp_size_t item_size_limbs;
    mp_size_t index_size_limbs;

    /* elliptic curve parameters a, b will 
    be represented in Montgomery form. */    
    mp_limb_t* curve_aR;
    mp_limb_t* curve_bR;
    mp_limb_t* curve_p;
    mp_limb_t* curve_P;
    mp_limb_t* G_order;

    /* sorry for all the 3 stars pointers... */
    mp_limb_t***      thread_tortoise_item_caches;  // n threads, each thread has m cache values.
    mp_limb_t***      thread_hare_item_caches;      // n threads, each thread has m cache values.
    mp_limb_t***      thread_tortoise_index_caches; // n threads, each thread has m cache indices.
    mp_limb_t***      thread_hare_index_caches;     // n threads, each thread has m cache indices.
    unsigned long***  thread_read_counters;         // n threads, each thread has n read counters (using n-1) to m write counters.
    unsigned long**   thread_write_counters;        // n threads, each thread has m write counters of m caches.

    mp_limb_t** random_tG_add_skG;                  // r random points of t*G + s*kG (using X, Y coordinate in Montgomery Form.)
    mp_limb_t** random_ts;                          // r random multipliers t, s

    /* other temporary stuffs for each thread */
    ecc_ptemp* thread_ecc_ptemps;
    mp_limb_t** thread_tmp_cache_reads;

    /* result storage */
    mp_limb_t** thread_result_tortoise_items;
    mp_limb_t** thread_result_hare_items;
    mp_limb_t** thread_result_tortoise_indices;
    mp_limb_t** thread_result_hare_indices;

    /* result indicators */
    int* founds;
    int overall_found;

} __dlog_obj_struct;

typedef __dlog_obj_struct dlog_obj[1];

int dlog_validate_input(
    ecc curve,
    eccpt G, eccpt kG,
    mpz_t G_mult_order,

    unsigned int n_threads,
    unsigned int n_caches,
    unsigned int n_randindices
);

int dlog_fast_solve_if_possible(
    ecc curve,
    mpz_t k, 
    eccpt G, eccpt kG,
    mpz_t G_mult_order
);

void dlog_init_dlog_obj(
    dlog_obj obj,

    ecc curve,
    eccpt G, eccpt kG,
    mpz_t G_mult_order,

    unsigned int n_threads,
    unsigned int n_caches,
    unsigned int n_randindices
);

void dlog_fill_dlog_obj(
    dlog_obj obj,

    ecc curve,
    eccpt G, eccpt kG,
    mpz_t G_mult_order,

    unsigned int n_threads,
    unsigned int n_caches,
    unsigned int n_randindices
);

void dlog_free_dlog_obj(
    dlog_obj obj
);

int dlog2(
    ecc curve, 
    mpz_t k, 
    eccpt G, eccpt kG, 
    mpz_t G_mult_order, 

    unsigned int n_threads,
    unsigned int n_caches,
    unsigned int n_randindices
);

#endif