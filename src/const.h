#ifndef CONST_H
#define CONST_H

#define AUTHOR                   "mistsuu"
#define BASESIZE_PRIME_CHECKER   40
#define UNDEFINED_COORDINATE     -1

#define DLOG_SUCCESS             0
#define DLOG_NOT_FOUND_DLOG      1
#define DLOG_MOVE_TO_NEXT_STEP   4
#define DLOG_POINT_NOT_ON_CURVE  5
#define DLOG_FAULTY_POINT_ORDER  6
#define DLOG_BAD_CONFIG          7
#define DLOG_BAD_COLLISION       8

#define DEFAULT_NUM_THREADS      4
#define DEFAULT_NUM_CACHE_ITEMS  4
#define DEFAULT_NUM_RAND_ITEMS   20

#define SHARED_LIB_NAME          "libdlogefp"

#endif