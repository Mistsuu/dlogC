#ifndef QUE2_H
#define QUE2_H

#include <gmp.h>
#include <stdlib.h>

typedef struct que2_struct
{
    signed long* data;
    signed long size;
    signed long max_size;
} __que2_struct;

typedef __que2_struct que2[1];

void que2_init(que2 queue);
void que2_push(que2 queue, signed long val1, signed long val2);
int que2_pop(que2 queue, signed long* val1, signed long* val2);
void que2_free(que2 queue);

#endif