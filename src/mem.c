#include "mem.h"

void* malloc_exit_when_null(size_t size)
{
    void* new_ptr = malloc(size);
    if (!new_ptr) {
        printf("[error] Error! Allocating memory gone wrong! Exiting...\n");
        exit(-1);
    }
    return new_ptr;
}