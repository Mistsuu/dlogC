#ifndef OLDDLOG_H
#define OLDDLOG_H

#include <gmp.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>

#ifdef DLOG_VERBOSE
    #include <sys/time.h>
#endif

#include "ecc.h"
#include "ecc_x.h"

size_t dlog_calc_mem(
    mpz_t n,
    size_t* index_size_limbs, size_t* index_size_bytes,
    size_t* item_size_limbs, size_t* item_size_bytes,
    size_t* n_partitions,
    
    mpz_t upper_k,
    size_t mem_limit,
    mpz_t curve_p
);

int dlog_alloc_buffer(
    char** lbuffer,
    char** rbuffer,

    size_t n_size_t, 
    size_t index_size_bytes, size_t item_size_bytes
);

typedef struct {
    char* _buffer;  

    size_t n_size_t; 
    size_t index_size_bytes; size_t index_size_limbs;
    size_t item_size_bytes; size_t item_size_limbs;

    mp_limb_t* _0x; mp_limb_t* _0z;   // must be modified in-use
    mp_limb_t* _1x; mp_limb_t* _1z;   // must be modified in-use
    mp_limb_t* dGx; mp_limb_t* dGz;

    mp_limb_t* i;
    int is_inc_i;

    mp_limb_t* curve_a;
    mp_limb_t* curve_b;
    mp_limb_t* curve_p;
} __args_thread__dlog_fill_buffer;

void* __thread__dlog_fill_buffer(
    void* vargs
);

void dlog_fill_buffer_l(
    char* lbuffer, 
    ecc curve, eccpt G, eccpt kG, 

    mpz_t n, size_t n_size_t, 
    size_t index_size_bytes, size_t index_size_limbs,
    size_t item_size_bytes, size_t item_size_limbs,
    
    unsigned int n_threads
);

void dlog_fill_buffer_r(
    char* rbuffer, 
    ecc curve, eccpt G, eccpt kG, 

    mpz_t n, size_t n_size_t, 
    size_t index_size_bytes, size_t index_size_limbs,
    size_t item_size_bytes, size_t item_size_limbs,
    
    unsigned int n_threads
);

void dlog_sort_buffer(
    char* _buffer,

    size_t n_size_t,
    size_t index_size_bytes,
    size_t item_size_bytes,

    unsigned int n_threads
);

typedef struct {
    mp_limb_t* exp_l_limbs;     // Must be index_size_limbs + 1 limbs allocated.
    mp_limb_t* exp_r_limbs;     // Must be index_size_limbs + 1 limbs allocated.

    char* lbuffer;
    char* rbuffer;

    size_t n_size_t_l;
    size_t n_size_t_r;

    size_t index_size_limbs; 
    size_t index_size_bytes;
    size_t item_size_bytes;
} __args_thread__dlog_search_buffer;

void* __thread__dlog_search_buffer(
    void* vargs
);

int dlog_search_buffer(
    mpz_t exp_l,
    mpz_t exp_r,

    char* lbuffer,
    char* rbuffer,
    
    size_t n_size_t, 
    size_t index_size_limbs, size_t index_size_bytes, 
    size_t item_size_bytes,

    unsigned int n_threads
);

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
);

int dlog(
    ecc curve, 
    mpz_t k, 
    eccpt G, eccpt kG, 
    mpz_t upper_k, 

    unsigned int n_threads,
    size_t mem_limit
);

#endif