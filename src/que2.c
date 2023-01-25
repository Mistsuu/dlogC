#include "que2.h"

void que2_init(que2* queue)
{
    queue->max_size = 4 * sizeof(size_t) * 8;
    queue->data = (size_t*) malloc(queue->max_size * sizeof(size_t));
    queue->size = 0;
}

void que2_push(que2* queue, size_t val1, size_t val2)
{
    if (queue->max_size > queue->size + 2) {
        queue->max_size += 4 * sizeof(size_t) * 8;
        queue->data = (size_t*) realloc(queue->data, queue->max_size * sizeof(size_t));
    }
    queue->data[queue->size++] = val1;
    queue->data[queue->size++] = val2;
}

int que2_pop(que2* queue, size_t* val1, size_t* val2) 
{
    if (queue->size < 2)
        return 0;
    (*val2) = queue->data[queue->size--];
    (*val1) = queue->data[queue->size--];
    return 1;
}

void que2_free(que2* queue)
{
    free(queue->data);
}