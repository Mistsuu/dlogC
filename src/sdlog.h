#ifndef SDLOG_H
#define SDLOG_H

#include <gmp.h>
#include <stdint.h>
#include <string.h>
#include <signal.h>

void user_interrupt_handler(
    int signum
);

int sdlog(
    // Field parameter
    char* str_p,

    // To be modified
    char** pstr_k,

    char* str_G,
    char* str_kG,
    char* str_n,
    
    // Configs
    unsigned long n_threads,
    unsigned long alpha,
    unsigned long n_rand_items
);

void sdlog_free(
    char* str_k
);

#endif