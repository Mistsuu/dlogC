#include "ecc_proj.h"
#include "ex_mpn.h"
#include "mem.h"

void ecc_init_ptemp(ecc_ptemp T, mp_size_t n)
{
    T[0]  = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
    T[1]  = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
    T[2]  = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
    T[3]  = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
    T[4]  = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
    T[5]  = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
    T[6]  = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
    T[7]  = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
    T[8]  = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
    T[9]  = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
    T[10] = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
    T[11] = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * (6*n));
}

void ecc_free_ptemp(ecc_ptemp T)
{
    free(T[0]);
    free(T[1]);
    free(T[2]);
    free(T[3]);
    free(T[4]);
    free(T[5]);
    free(T[6]);
}

// -------------------------------------------------------------------------------

void ecc_padd(
    mp_limb_t* Rx, mp_limb_t* Ry, mp_limb_t* Rz,       // Rx, Ry, Rz must have n limbs allocated
    mp_limb_t* Px, mp_limb_t* Py, mp_limb_t* Pz,       // Px, Ry, Pz must have n limbs allocated
    mp_limb_t* Qx, mp_limb_t* Qy,                      // Qx, Qy     must have n limbs allocated

    mp_limb_t* curve_p,                                // curve_p must have n limbs allocated
    mp_limb_t* curve_P,                                // curve_P must have n limbs allocated
    mp_size_t n,                                       // number of limbs in curve->p
    
    ecc_ptemp T                                        // temporary variables, allocated with ecc_init_ptemp(T, n).
)
{
    // u = Y2*Z1
    mpn_montgomery_mulmod_n(T[0], Qy, Pz, curve_p, curve_P, n, T[11]);
    // u = u-Y1
    mpn_montgomery_submod_n(T[0], T[0], Py, curve_p, n);
    // uu = u**2
    mpn_montgomery_sqrmod_n(T[1], T[0], curve_p, curve_P, n, T[11]);
    // v = X2*Z1
    mpn_montgomery_mulmod_n(T[2], Qx, Pz, curve_p, curve_P, n, T[11]);
    // v = v-X1
    mpn_montgomery_submod_n(T[2], T[2], Px, curve_p, n);
    // vv = v**2
    mpn_montgomery_sqrmod_n(T[3], T[2], curve_p, curve_P, n, T[11]);
    // vvv = v*vv
    mpn_montgomery_mulmod_n(T[4], T[2], T[3], curve_p, curve_P, n, T[11]);
    // R = vv*X1
    mpn_montgomery_mulmod_n(T[5], T[3], Px, curve_p, curve_P, n, T[11]);
    // A = uu*Z1
    mpn_montgomery_mulmod_n(T[6], T[1], Pz, curve_p, curve_P, n, T[11]);
    // A = A-vvv
    mpn_montgomery_submod_n(T[6], T[6], T[4], curve_p, n);
    // A = A-R
    mpn_montgomery_submod_n(T[6], T[6], T[5], curve_p, n);
    // A = A-R
    mpn_montgomery_submod_n(T[6], T[6], T[5], curve_p, n);
    // X3 = v*A
    mpn_montgomery_mulmod_n(Rx, T[2], T[6], curve_p, curve_P, n, T[11]);
    // t5 = R-A
    mpn_montgomery_submod_n(T[7], T[5], T[6], curve_p, n);
    // t6 = vvv*Y1
    mpn_montgomery_mulmod_n(T[8], T[4], Py, curve_p, curve_P, n, T[11]);
    // Y3 = u*t5
    mpn_montgomery_mulmod_n(Ry, T[0], T[7], curve_p, curve_P, n, T[11]);
    // Y3 = Y3-t6
    mpn_montgomery_submod_n(Ry, Ry, T[8], curve_p, n);
    // Z3 = vvv*Z1
    mpn_montgomery_mulmod_n(Rz, T[4], Pz, curve_p, curve_P, n, T[11]);
}

void ecc_pdbl(
    mp_limb_t* Rx, mp_limb_t* Ry, mp_limb_t* Rz,       // Rx, Ry, Rz must have n limbs allocated
    mp_limb_t* Px, mp_limb_t* Py, mp_limb_t* Pz,       // Px, Ry, Pz must have n limbs allocated

    mp_limb_t* curve_aR,                               // curve_aR must have n limbs allocated
    mp_limb_t* curve_p,                                // curve_p must have n limbs allocated
    mp_limb_t* curve_P,                                // curve_P must have n limbs allocated
    mp_size_t n,                                       // number of limbs in curve->p
    
    ecc_ptemp T                                        // temporary variables, allocated with ecc_init_ptemp(T, n).
)
{
    // XX = X1**2
    mpn_montgomery_sqrmod_n(T[0], Px, curve_p, curve_P, n, T[11]);
    // ZZ = Z1**2
    mpn_montgomery_sqrmod_n(T[1], Pz, curve_p, curve_P, n, T[11]);
    // t0 = XX+XX
    mpn_montgomery_lshift1mod_n(T[2], T[0], curve_p, n);
    // t0 = t0+XX
    mpn_montgomery_addmod_n(T[2], T[2], T[0], curve_p, n);
    // w = a*ZZ
    mpn_montgomery_mulmod_n(T[3], curve_aR, T[1], curve_p, curve_P, n, T[11]);
    // w = w+t0
    mpn_montgomery_addmod_n(T[3], T[3], T[2], curve_p, n);
    // s = Y1*Z1
    mpn_montgomery_mulmod_n(T[4], Py, Pz, curve_p, curve_P, n, T[11]);
    // s = s+s
    mpn_montgomery_lshift1mod_n(T[4], T[4], curve_p, n);
    // ss = s**2
    mpn_montgomery_sqrmod_n(T[5], T[4], curve_p, curve_P, n, T[11]);
    // Z3 = s*ss
    mpn_montgomery_mulmod_n(Rz, T[4], T[5], curve_p, curve_P, n, T[11]);
    // R = Y1*s
    mpn_montgomery_mulmod_n(T[6], Py, T[4], curve_p, curve_P, n, T[11]);
    // RR = R**2
    mpn_montgomery_sqrmod_n(T[7], T[6], curve_p, curve_P, n, T[11]);
    // t3 = X1+R
    mpn_montgomery_addmod_n(T[8], Px, T[6], curve_p, n);
    // B = t3**2
    mpn_montgomery_sqrmod_n(T[9], T[8], curve_p, curve_P, n, T[11]);
    // B = B-XX
    mpn_montgomery_submod_n(T[9], T[9], T[0], curve_p, n);
    // B = B-RR
    mpn_montgomery_submod_n(T[9], T[9], T[7], curve_p, n);
    // h = w**2
    mpn_montgomery_sqrmod_n(T[10], T[3], curve_p, curve_P, n, T[11]);
    // h = h-B
    mpn_montgomery_submod_n(T[10], T[10], T[9], curve_p, n);
    // h = h-B
    mpn_montgomery_submod_n(T[10], T[10], T[9], curve_p, n);
    // X3 = h*s
    mpn_montgomery_mulmod_n(Rx, T[10], T[4], curve_p, curve_P, n, T[11]);
    // Y3 = B-h
    mpn_montgomery_submod_n(Ry, T[9], T[10], curve_p, n);
    // Y3 = w*Y3
    mpn_montgomery_mulmod_n(Ry, T[3], Ry, curve_p, curve_P, n, T[11]);
    // Y3 = Y3-RR
    mpn_montgomery_submod_n(Ry, Ry, T[7], curve_p, n);
    // Y3 = Y3-RR
    mpn_montgomery_submod_n(Ry, Ry, T[7], curve_p, n);
}