#include "que2.h"
#include "mem.h"

void que2_init(que2 queue)
{
    queue->max_size = 4 * sizeof(signed long) * 8;
    queue->data = (signed long*) malloc_exit_when_null(queue->max_size * sizeof(signed long));
    queue->size = 0;
}

void que2_push(que2 queue, signed long val1, signed long val2)
{
    if (queue->max_size > queue->size + 2) {
        queue->max_size += 4 * sizeof(signed long) * 8;
        queue->data = (signed long*) realloc(queue->data, queue->max_size * sizeof(signed long));
    }
    queue->data[queue->size++] = val1;
    queue->data[queue->size++] = val2;
}

int que2_pop(que2 queue, signed long* val1, signed long* val2) 
{
    if (queue->size < 2)
        return 0;
    (*val2) = queue->data[--queue->size];
    (*val1) = queue->data[--queue->size];
    return 1;
}

void que2_free(que2 queue)
{
    free(queue->data);
}