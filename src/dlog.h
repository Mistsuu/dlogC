#ifndef DLOG_H
#define DLOG_H

#include <gmp.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <xxhash.h>

#ifdef DLOG_VERBOSE
    #include <sys/time.h>
#endif

#include "ecc.h"
#include "ecc_proj.h"

typedef struct dlog_obj_struct
{   
    unsigned int n_threads;
    size_t mem_limit;
    unsigned int n_rand_items;
    mp_size_t item_size_limbs;
    mp_size_t index_size_limbs;

    mp_size_t hash_item_size_limbs;
    size_t n_hash_items;

    /* elliptic curve parameters a, b will 
    be represented in Montgomery form. */    
    mp_limb_t* curve_aR;
    mp_limb_t* curve_bR;
    mp_limb_t* curve_p;
    mp_limb_t* curve_P;
    mp_limb_t* G_order;

    /* okay, this version doesn't have any 3 stars pointers :> */
    mp_limb_t** thread_tortoise_X_items;          // n threads, each thread has 1 value.
    mp_limb_t** thread_tortoise_ts_indices;       // n threads, each thread has 1 value.
    mp_limb_t** thread_hare_XYZ_items;            // n threads, each thread has 1 value.
    mp_limb_t** thread_hare_X_items;              // n threads, each thread has 1 value.
    mp_limb_t** thread_hare_ts_index;             // n threads, each thread has 1 value.

    /* results points are hashed and t,s indices are put into the collision area... */
    mp_limb_t*  ts_index_hashstores;              // place to store t, s indices.

    /* random points used in Teske's method to quickly generate psuedo-random points */
    mp_limb_t** random_tG_add_skG;                // r random points of t*G + s*kG (using X, Y coordinate in Montgomery Form.)
    mp_limb_t** random_ts;                        // r random multipliers t, s

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
    ecc curve,
    eccpt G, eccpt kG,
    mpz_t G_mult_order,

    unsigned int n_threads,
    size_t mem_limit,
    unsigned int n_rand_items
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
    size_t mem_limit,
    unsigned int n_rand_items
);

void dlog_fill_dlog_obj(
    dlog_obj obj,

    ecc curve,
    eccpt G, eccpt kG,
    mpz_t G_mult_order,

    unsigned int n_threads,
    size_t mem_limit,
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
    ecc curve,
    mpz_t k,
    eccpt G, eccpt kG,
    mpz_t G_mult_order,

    dlog_obj obj
);

int dlog(
    ecc curve, 
    mpz_t k, 
    eccpt G, eccpt kG, 
    mpz_t G_mult_order, 

    unsigned int n_threads,
    size_t mem_limit,
    unsigned int n_rand_items
);

#endif