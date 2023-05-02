#include "num.h"
#include "mem.h"

void mpn_printf(mp_limb_t* mpn_in, mp_size_t mpn_len)
{
    mpz_t mpz_in;
    mpz_in->_mp_alloc = mpn_len;
    mpz_in->_mp_size  = mpn_len;
    mpz_in->_mp_d     = mpn_in;
    mpz_out_str(stdout, 10, mpz_in);
}

mp_limb_t* mpz_limbs_init_cpy(mpz_t x, size_t n)
{
    const mp_limb_t* x_limbs = mpz_limbs_read(x);
    size_t           x_size  = mpz_size(x);

    mp_limb_t* x_copied_limbs = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
    memset(x_copied_limbs, 0, sizeof(mp_limb_t) * n);
    memcpy(x_copied_limbs, x_limbs, sizeof(mp_limb_t) * (x_size > n ? n : x_size));
    return x_copied_limbs;
}

mp_limb_t* mpz_limbs_init_zero(size_t n)
{
    mp_limb_t* x_limbs = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
    memset(x_limbs, 0, sizeof(mp_limb_t) * n);
    return x_limbs;
}

void mpn2bytes(unsigned char *str, mp_size_t len, const mp_limb_t *s1p, mp_size_t s1n)
{
    // If I don't do this, 
    // mpn_get_str() will return some crazy stuffs :'<
    if (mpn_zero_p(s1p, s1n)) {
        memset(str, 0, len);
        return;
    }

    mp_size_t actual_len = mpn_get_str(str, 256, (mp_limb_t *)s1p, s1n);
    if (actual_len > len) { // todo: god i hope i could delete this shit
        printf("[error] omg plz, whyyyyy mpn_get_str cannot fit in mpn2bytes: actual_len=%ld, len=%ld???????\n", actual_len, len);
        printf("[error] s1p: \n");

        for (mp_size_t i = 0; i < s1n; ++i) {
            printf("       0x%016lx\n", s1p[i]);
        }
        exit(-1);
    }
    
    memmove(str + len - actual_len, str, actual_len);
    memset(str, 0, len - actual_len);
}