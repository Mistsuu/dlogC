#ifndef DLOG_H
#define DLOG_H

#include <gmp.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>

#ifdef DLOG_VERBOSE
    #include <time.h>
#endif

#include "ecc.h"
#include "ecc_x.h"

size_t dlog_init_buffer(
    char** lbuffer,
    char** rbuffer,
    mpz_t n, size_t index_size_bytes, size_t item_size_bytes
);

typedef struct {
    char* lbuffer;  
    char* rbuffer;

    size_t n_size_t; 
    size_t index_size_bytes; size_t index_size_limbs;
    size_t item_size_bytes; size_t item_size_limbs;

    mp_limb_t* L0x; mp_limb_t* L0z;   // must be modified in-use
    mp_limb_t* L1x; mp_limb_t* L1z;   // must be modified in-use
    mp_limb_t* R0x; mp_limb_t* R0z;   // must be modified in-use
    mp_limb_t* R1x; mp_limb_t* R1z;   // must be modified in-use
    mp_limb_t* Gx; mp_limb_t* Gz;
    mp_limb_t* nGx; mp_limb_t* nGz;

    mp_limb_t* i_l; mp_limb_t* i_r;

    mp_limb_t* curve_a;
    mp_limb_t* curve_b;
    mp_limb_t* curve_p;
} __args_thread__dlog_fill_buffer;

void* __thread__dlog_fill_buffer(
    void* vargs
);

void dlog_fill_buffer(
    char* lbuffer,
    char* rbuffer,
    ecc curve, eccpt G, eccpt kG, 

    mpz_t n, size_t n_size_t, 
    size_t index_size_bytes, size_t index_size_limbs,
    size_t item_size_bytes, size_t item_size_limbs,
    
    unsigned int n_threads
);

void dlog_sort_buffer(
    char* lbuffer,
    char* rbuffer,

    size_t n_size_t,
    size_t index_size_bytes,
    size_t item_size_bytes
);

int dlog_search_buffer(
    mpz_t exp_l,
    mpz_t exp_r,

    char* lbuffer,
    char* rbuffer,
    
    size_t n_size_t, 
    size_t index_size_limbs, size_t index_size_bytes, 
    size_t item_size_bytes
);

int dlog(ecc curve, mpz_t k, eccpt G, eccpt kG, mpz_t upper_k, unsigned int n_threads);

#endif