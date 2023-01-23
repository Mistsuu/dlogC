#include "gmp.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char * argv[]) {
    printf("%s\n", gmp_version);
    printf("%d\n", mp_bits_per_limb);
}