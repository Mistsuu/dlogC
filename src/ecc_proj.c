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
    T[7]  = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * (6*n));
    T[8]  = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
    T[9]  = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
    T[10] = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * n);
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
    free(T[7]);
}

// -------------------------------------------------------------------------------

void ecc_padd(
    mp_limb_t* Rx, mp_limb_t* Ry, mp_limb_t* Rz,       // Rx, Ry, Rz must have n limbs allocated
    mp_limb_t* Px, mp_limb_t* Py, mp_limb_t* Pz,       // Px, Ry, Pz must have n limbs allocated
    mp_limb_t* Qx, mp_limb_t* Qy, mp_limb_t* Qz,       // Qx, Qy, Qz must have n limbs allocated

    mp_limb_t* curve_a,                                // curve_a must have n limbs allocated
    mp_limb_t* curve_b3,                               // curve_b3 must have n limbs allocated
    mp_limb_t* curve_p,                                // curve_p must have n limbs allocated
    mp_limb_t* curve_P,                                // curve_P must have n limbs allocated
    mp_size_t n,                                       // number of limbs in curve->p
    
    ecc_ptemp T                                        // temporary variables, allocated with ecc_init_ptemp(T, n).
)
{
    // t0 = X1*X2
    mpn_montgomery_mulmod_n(T[0], Px, Qx, curve_p, curve_P, n, T[7]);
    // t1 = Y1*Y2
    mpn_montgomery_mulmod_n(T[1], Py, Qy, curve_p, curve_P, n, T[7]);
    // t2 = Z1*Z2
    mpn_montgomery_mulmod_n(T[2], Pz, Qz, curve_p, curve_P, n, T[7]);
    // t3 = X1+Y1
    mpn_montgomery_addmod_n(T[3], Px, Py, curve_p, n);
    // t4 = X2+Y2
    mpn_montgomery_addmod_n(T[4], Qx, Qy, curve_p, n);
    // t3 = t3*t4
    mpn_montgomery_mulmod_n(T[6], T[3], T[4], curve_p, curve_P, n, T[7]);
    mpn_copyd(T[3], T[6], n);
    // t4 = t0+t1
    mpn_montgomery_addmod_n(T[4], T[0], T[1], curve_p, n);
    // t3 = t3-t4
    mpn_montgomery_submod_n(T[3], T[3], T[4], curve_p, n);
    // t4 = X1+Z1
    mpn_montgomery_addmod_n(T[4], Px, Pz, curve_p, n);
    // t5 = X2+Z2
    mpn_montgomery_addmod_n(T[5], Qx, Qz, curve_p, n);
    // t4 = t4*t5
    mpn_montgomery_mulmod_n(T[6], T[4], T[5], curve_p, curve_P, n, T[7]);
    mpn_copyd(T[4], T[6], n);
    // t5 = t0+t2
    mpn_montgomery_addmod_n(T[5], T[0], T[2], curve_p, n);
    // t4 = t4-t5
    mpn_montgomery_submod_n(T[4], T[4], T[5], curve_p, n);
    // t5 = Y1+Z1
    mpn_montgomery_addmod_n(T[5], Py, Pz, curve_p, n);
    // X3 = Y2+Z2
    mpn_montgomery_addmod_n(Rx, Qy, Qz, curve_p, n);
    // t5 = t5*X3
    mpn_montgomery_mulmod_n(T[6], T[5], Rx, curve_p, curve_P, n, T[7]);
    mpn_copyd(T[5], T[6], n);
    // X3 = t1+t2
    mpn_montgomery_addmod_n(Rx, T[1], T[2], curve_p, n);
    // t5 = t5-X3
    mpn_montgomery_submod_n(T[5], T[5], Rx, curve_p, n);
    // Z3 = a*t4
    mpn_montgomery_mulmod_n(Rz, curve_a, T[4], curve_p, curve_P, n, T[7]);
    // X3 = b3*t2
    mpn_montgomery_mulmod_n(Rx, curve_b3, T[2], curve_p, curve_P, n, T[7]);
    // Z3 = X3+Z3
    mpn_montgomery_addmod_n(Rz, Rx, Rz, curve_p, n);
    // X3 = t1-Z3
    mpn_montgomery_submod_n(Rx, T[1], Rz, curve_p, n);
    // Z3 = t1+Z3
    mpn_montgomery_addmod_n(Rz, T[1], Rz, curve_p, n);
    // Y3 = X3*Z3
    mpn_montgomery_mulmod_n(Ry, Rx, Rz, curve_p, curve_P, n, T[7]);
    // t1 = t0+t0
    mpn_montgomery_lshift1mod_n(T[1], T[0], curve_p, n);
    // t1 = t1+t0
    mpn_montgomery_addmod_n(T[1], T[1], T[0], curve_p, n);
    // t2 = a*t2
    mpn_montgomery_mulmod_n(T[2], curve_a, T[2], curve_p, curve_P, n, T[7]);
    // t4 = b3*t4
    mpn_montgomery_mulmod_n(T[4], curve_b3, T[4], curve_p, curve_P, n, T[7]);
    // t1 = t1+t2
    mpn_montgomery_addmod_n(T[1], T[1], T[2], curve_p, n);
    // t2 = t0-t2
    mpn_montgomery_submod_n(T[2], T[0], T[2], curve_p, n);
    // t2 = a*t2
    mpn_montgomery_mulmod_n(T[2], curve_a, T[2], curve_p, curve_P, n, T[7]);
    // t4 = t4+t2
    mpn_montgomery_addmod_n(T[4], T[4], T[2], curve_p, n);
    // t0 = t1*t4
    mpn_montgomery_mulmod_n(T[0], T[1], T[4], curve_p, curve_P, n, T[7]);
    // Y3 = Y3+t0
    mpn_montgomery_addmod_n(Ry, Ry, T[0], curve_p, n);
    // t0 = t5*t4
    mpn_montgomery_mulmod_n(T[0], T[5], T[4], curve_p, curve_P, n, T[7]);
    // X3 = t3*X3
    mpn_montgomery_mulmod_n(Rx, T[3], Rx, curve_p, curve_P, n, T[7]);
    // X3 = X3-t0
    mpn_montgomery_submod_n(Rx, Rx, T[0], curve_p, n);
    // t0 = t3*t1
    mpn_montgomery_mulmod_n(T[0], T[3], T[1], curve_p, curve_P, n, T[7]);
    // Z3 = t5*Z3
    mpn_montgomery_mulmod_n(Rz, T[5], Rz, curve_p, curve_P, n, T[7]);
    // Z3 = Z3+t0
    mpn_montgomery_addmod_n(Rz, Rz, T[0], curve_p, n);
}

void ecc_pdbl(
    mp_limb_t* Rx, mp_limb_t* Ry, mp_limb_t* Rz,       // Rx, Ry, Rz must have n limbs allocated
    mp_limb_t* Px, mp_limb_t* Py, mp_limb_t* Pz,       // Px, Ry, Pz must have n limbs allocated

    mp_limb_t* curve_a,                                // curve_a must have n limbs allocated
    mp_limb_t* curve_b3,                               // curve_b3 must have n limbs allocated
    mp_limb_t* curve_p,                                // curve_p must have n limbs allocated
    mp_limb_t* curve_P,                                // curve_P must have n limbs allocated
    mp_size_t n,                                       // number of limbs in curve->p
    
    ecc_ptemp T                                        // temporary variables, allocated with ecc_init_ptemp(T, n).
)
{
    // t0 = X1**2
    mpn_montgomery_sqrmod_n(T[0], Px, curve_p, curve_P, n, T[7]);
    // t1 = Y1**2
    mpn_montgomery_sqrmod_n(T[1], Py, curve_p, curve_P, n, T[7]);
    // t2 = Z1**2
    mpn_montgomery_sqrmod_n(T[2], Pz, curve_p, curve_P, n, T[7]);
    // t3 = X1*Y1
    mpn_montgomery_mulmod_n(T[3], Px, Py, curve_p, curve_P, n, T[7]);
    // t3 = t3+t3
    mpn_montgomery_lshift1mod_n(T[3], T[3], curve_p, n);
    // Z3 = X1*Z1
    mpn_montgomery_mulmod_n(Rz, Px, Pz, curve_p, curve_P, n, T[7]);
    // Z3 = Z3+Z3
    mpn_montgomery_lshift1mod_n(Rz, Rz, curve_p, n);
    // X3 = a*Z3
    mpn_montgomery_mulmod_n(Rx, curve_a, Rz, curve_p, curve_P, n, T[7]);
    // Y3 = b3*t2
    mpn_montgomery_mulmod_n(Ry, curve_b3, T[2], curve_p, curve_P, n, T[7]);
    // Y3 = X3+Y3
    mpn_montgomery_addmod_n(Ry, Rx, Ry, curve_p, n);
    // X3 = t1-Y3
    mpn_montgomery_submod_n(Rx, T[1], Ry, curve_p, n);
    // Y3 = t1+Y3
    mpn_montgomery_addmod_n(Ry, T[1], Ry, curve_p, n);
    // Y3 = X3*Y3
    mpn_montgomery_mulmod_n(Ry, Rx, Ry, curve_p, curve_P, n, T[7]);
    // X3 = t3*X3
    mpn_montgomery_mulmod_n(Rx, T[3], Rx, curve_p, curve_P, n, T[7]);
    // Z3 = b3*Z3
    mpn_montgomery_mulmod_n(Rz, curve_b3, Rz, curve_p, curve_P, n, T[7]);
    // t2 = a*t2
    mpn_montgomery_mulmod_n(T[2], curve_a, T[2], curve_p, curve_P, n, T[7]);
    // t3 = t0-t2
    mpn_montgomery_submod_n(T[3], T[0], T[2], curve_p, n);
    // t3 = a*t3
    mpn_montgomery_mulmod_n(T[3], curve_a, T[3], curve_p, curve_P, n, T[7]);
    // t3 = t3+Z3
    mpn_montgomery_addmod_n(T[3], T[3], Rz, curve_p, n);
    // Z3 = t0+t0
    mpn_montgomery_lshift1mod_n(Rz, T[0], curve_p, n);
    // t0 = Z3+t0
    mpn_montgomery_addmod_n(T[0], Rz, T[0], curve_p, n);
    // t0 = t0+t2
    mpn_montgomery_addmod_n(T[0], T[0], T[2], curve_p, n);
    // t0 = t0*t3
    mpn_montgomery_mulmod_n(T[6], T[0], T[3], curve_p, curve_P, n, T[7]);
    mpn_copyd(T[0], T[6], n);
    // Y3 = Y3+t0
    mpn_montgomery_addmod_n(Ry, Ry, T[0], curve_p, n);
    // t2 = Y1*Z1
    mpn_montgomery_mulmod_n(T[2], Py, Pz, curve_p, curve_P, n, T[7]);
    // t2 = t2+t2
    mpn_montgomery_lshift1mod_n(T[2], T[2], curve_p, n);
    // t0 = t2*t3
    mpn_montgomery_mulmod_n(T[0], T[2], T[3], curve_p, curve_P, n, T[7]);
    // X3 = X3-t0
    mpn_montgomery_submod_n(Rx, Rx, T[0], curve_p, n);
    // Z3 = t2*t1
    mpn_montgomery_mulmod_n(Rz, T[2], T[1], curve_p, curve_P, n, T[7]);
    // Z3 = Z3+Z3
    mpn_montgomery_lshift1mod_n(Rz, Rz, curve_p, n);
    // Z3 = Z3+Z3
    mpn_montgomery_lshift1mod_n(Rz, Rz, curve_p, n);
}

void ecc_pmul(
    mp_limb_t* Rx, mp_limb_t* Ry, mp_limb_t* Rz,       // Rx, Ry, Rz must have n limbs allocated
    mp_limb_t* Px, mp_limb_t* Py, mp_limb_t* Pz,       // Px, Ry, Pz must have n limbs allocated

    mp_limb_t* k,                                      // k must have n limbs allocated

    mp_limb_t* curve_a,                                // curve_a must have n limbs allocated
    mp_limb_t* curve_b3,                               // curve_b3 must have n limbs allocated
    mp_limb_t* curve_p,                                // curve_p must have n limbs allocated
    mp_limb_t* curve_P,                                // curve_p must have n limbs allocated
    mp_size_t n,                                       // number of limbs in curve->p
    
    ecc_ptemp T                                        // temporary variables, allocated with ecc_init_ptemp(T, n).
)
{
    int i = n-1;
    int j = 0;
    mp_limb_t k_limb;

    /* ---------- Find the first bit of k != 0 ---------- */
    while (i>=0 && !k[i])
        i--;

    if (i < 0) {
        mpn_zero(Rx, n);
        mpn_zero(Ry, n);
        mpn_zero(Rz, n);
        return;
    }

    k_limb = k[i];
    while (k_limb) {
        k_limb >>= 1;
        j++;
    }

    j -= 2;

    /* ---------- Initialize ---------- */
    // R0 = P
    mpn_copyi(Rx, Px, n);
    mpn_copyi(Ry, Py, n);
    mpn_copyi(Rz, Pz, n);

    // R1 = P*2
    mp_limb_t* Tx = T[8];
    mp_limb_t* Ty = T[9];
    mp_limb_t* Tz = T[10];
    ecc_pdbl(
        Tx, Ty, Tz,
        Px, Py, Pz,

        curve_a,
        curve_b3,
        curve_p,
        curve_P,
        n,

        T
    );

    /* ---------- Doing ladder ---------- */
    while (i>=0) {
        k_limb = k[i];
        while (j>=0) {
            if ((k_limb >> j) & 1) {
                // R0 <- R1 + R0
                ecc_padd(
                    Rx, Ry, Rz,
                    Tx, Ty, Tz,
                    Px, Py, Pz,

                    curve_a,
                    curve_b3,
                    curve_p,
                    curve_P,
                    n,

                    T
                );

                // R1 <- 2*R1
                ecc_pdbl(
                    Tx, Ty, Tz,
                    Tx, Ty, Tz,

                    curve_a,
                    curve_b3,
                    curve_p,
                    curve_P,
                    n,

                    T
                );
            }
            else {
                // R1 <- R1 + R0
                ecc_padd(
                    Tx, Ty, Tz,
                    Rx, Ry, Rz,
                    Px, Py, Pz,

                    curve_a,
                    curve_b3,
                    curve_p,
                    curve_P,
                    n,

                    T
                );

                // R0 <- 2*R0
                ecc_pdbl(
                    Rx, Ry, Rz,
                    Rx, Ry, Rz,

                    curve_a,
                    curve_b3,
                    curve_p,
                    curve_P,
                    n,

                    T
                );
            }
            j--;
        }
        i--;
        j = mp_bits_per_limb-1;
    }

    // R <- R0 (already done)
}