#ifndef SDLOG_H
#define SDLOG_H

#include <gmp.h>
#include <stdint.h>
#include <string.h>

int sdlog(
    char* str_p,
    char** pstr_k,
    char* str_G,
    char* str_kG,
    char* str_upper_k,
    unsigned int n_threads,
    size_t mem_limit
);

void sdlog_free(
    char* str_k
);

#endif