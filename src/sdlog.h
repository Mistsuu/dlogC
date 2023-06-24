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
    // Curve parameters
    char* str_curve_a,
    char* str_curve_b,
    char* str_curve_p,

    // To be overwritten
    char** pstr_k,

    char* str_Gx,
    char* str_Gy,
    char* str_Gz,
    char* str_kGx,
    char* str_kGy,
    char* str_kGz,
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