#ifndef QUE2_H
#define QUE2_H

#include <gmp.h>
#include <stdlib.h>

typedef struct que2_struct
{
    size_t* data;
    size_t size;
    size_t max_size;
} __que2_struct;

typedef __que2_struct que2[1];

void que2_init(que2 queue);
void que2_push(que2 queue, size_t val1, size_t val2);
int que2_pop(que2 queue, size_t* val1, size_t* val2);
void que2_free(que2 queue);

#endif