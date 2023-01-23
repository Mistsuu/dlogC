#include "num.h"

void mpn_printf(mp_limb_t* mpn_in, mp_size_t mpn_len)
{
    mpz_t mpz_in;
    mpz_in->_mp_alloc = mpn_len;
    mpz_in->_mp_size  = mpn_len;
    mpz_in->_mp_d     = mpn_in;
    mpz_out_str(stdout, 10, mpz_in);
}