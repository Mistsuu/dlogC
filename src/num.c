#include "num.h"

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

    mp_limb_t* x_copied_limbs = (mp_limb_t*) malloc(sizeof(mp_limb_t) * n);
    memset(x_copied_limbs, 0, sizeof(mp_limb_t) * n);
    memcpy(x_copied_limbs, x_limbs, sizeof(mp_limb_t) * (x_size > n ? n : x_size));
    return x_copied_limbs;
}