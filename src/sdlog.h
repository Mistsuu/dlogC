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
    char* str_upper_k,
    
    // Configs
    unsigned int n_threads,
    size_t mem_limit
);

void sdlog_free(
    char* str_k
);

#endif