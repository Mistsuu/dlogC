#ifndef DLOG_H
#define DLOG_H

#include <gmp.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>

#ifdef DLOG_VERBOSE
    #include <sys/time.h>
#endif

size_t dlog_calc_mem(
    mpz_t n,
    size_t* index_size_limbs, size_t* index_size_bytes,
    size_t* item_size_limbs, size_t* item_size_bytes,
    size_t* n_partitions,
    
    mpz_t upper_k,
    size_t mem_limit,
    mpz_t p
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

    mp_limb_t* dG;
    mp_limb_t* _0;
    mp_limb_t* p;

    mp_limb_t* i;
    int is_inc_i;
} __args_thread__dlog_fill_buffer;

void* __thread__dlog_fill_buffer(
    void* vargs
);

void dlog_fill_buffer_l(
    char* lbuffer, 
    mpz_t p, mpz_t G, mpz_t kG, 

    mpz_t n, size_t n_size_t, 
    size_t index_size_bytes, size_t index_size_limbs,
    size_t item_size_bytes, size_t item_size_limbs,
    
    unsigned int n_threads
);

void dlog_fill_buffer_r(
    char* rbuffer, 
    mpz_t p, mpz_t G, mpz_t kG,

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
    mpz_t p,
    mpz_t k, 
    mpz_t G, 
    mpz_t kG, 
    
    char* lbuffer, 
    char* rbuffer,
    
    mpz_t n, size_t n_size_t,
    size_t index_size_limbs, size_t index_size_bytes,
    size_t item_size_limbs, size_t item_size_bytes,

    unsigned int n_threads,
    unsigned int is_update_lbuffer
);

int dlog(
    mpz_t p, 
    mpz_t k, 
    mpz_t G, 
    mpz_t kG, 

    mpz_t upper_k, 

    unsigned int n_threads,
    size_t mem_limit
);

#endif