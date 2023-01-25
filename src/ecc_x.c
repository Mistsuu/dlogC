#include "ecc_x.h"
#include "num.h"
#include "mem.h"

void ecc_init_xtemp(ecc_xtemp T, mp_size_t n)
{
    T[0] = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * (5*n+2));
    T[1] = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * (2*n+1));
    T[2] = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * (2*n));
    T[3] = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * (4*n));
    T[4] = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * (6*n));
    T[5] = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * (5*n+2));
    T[6] = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * (5*n+1));
    T[7] = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * (n));
    T[8] = (mp_limb_t*) malloc_exit_when_null(sizeof(mp_limb_t) * (n));
}

void ecc_free_xtemp(ecc_xtemp T)
{
    free(T[0]);
    free(T[1]);
    free(T[2]);
    free(T[3]);
    free(T[4]);
    free(T[5]);
    free(T[6]);
    free(T[7]);
    free(T[8]);
}

// -------------------------------------------------------------------------------

void ecc_xadd(
    mp_limb_t* Rx, mp_limb_t* Rz,       // Rx, Rz must have n limbs allocated
    mp_limb_t* Px, mp_limb_t* Pz,       // Px, Pz must have n limbs allocated
    mp_limb_t* Qx, mp_limb_t* Qz,       // Qx, Qz must have n limbs allocated
    mp_limb_t* P_Qx, mp_limb_t* P_Qz,   // P_Qx, P_Qz must have n limbs allocated

    mp_limb_t* curve_a,                 // curve_a must have n limbs allocated
    mp_limb_t* curve_b,                 // curve_b must have n limbs allocated
    mp_limb_t* curve_p,                 // curve_p must have n limbs allocated
    mp_size_t n,                        // number of limbs in curve->p
    
    ecc_xtemp T                         // temporary variables, allocated with ecc_init_xtemp(T, n).
)
{
    // T0 = Px*Qx
    mpn_mul_n(T[0], Px, Qx, n);
    // T1 = Pz*Qz
    mpn_mul_n(T[1], Pz, Qz, n);
    // T2 = Px*Qz
    mpn_mul_n(T[2], Px, Qz, n);
    // T3 = Pz*Qx
    mpn_mul_n(T[3], Pz, Qx, n);
    // T4 = a*T1
    mpn_mul(T[4], T[1], 2*n, curve_a, n);
    // T0 = |T0-T4|
    if (mpn_sub(T[0], T[4], 3*n, T[0], 2*n))
        mpn_neg(T[0], T[0], 3*n);
    // T4 = T0^2
    mpn_sqr(T[4], T[0], 3*n);
    // T0 = b*T1
    mpn_mul(T[0], T[1], 2*n, curve_b, n);
    // T0 = 4*T0
    T[0][3*n] = mpn_lshift(T[0], T[0], 3*n, 2);
    // T1 = T2+T3
    T[1][2*n] = mpn_add_n(T[1], T[2], T[3], 2*n);
    // T5 = T0*T1
    mpn_mul(T[5], T[0], 3*n+1, T[1], 2*n+1);
    // T4 = (T4-T5) mod p
    if (mpn_sub(T[4], T[4], 6*n, T[5], 5*n+2) == 0) {
        // T4 >= T5: T4 = (T4-T5) mod p
        mpn_tdiv_qr(T[6], T[4], 0, T[4], 6*n, curve_p, n);
    } else {
        // T4 < T5: T4 = p - (T5-T4) mod p
        mpn_neg(T[4], T[4], 6*n);
        mpn_tdiv_qr(T[6], T[4], 0, T[4], 6*n, curve_p, n);
        mpn_sub_n(T[4], curve_p, T[4], n);
    }
    // Rx = P_Qz*T4
    mpn_mul_n(T[5], P_Qz, T[4], n);
    mpn_tdiv_qr(T[6], Rx, 0, T[5], 2*n, curve_p, n);
    // T2 = |T2-T3|
    if (mpn_sub_n(T[2], T[2], T[3], 2*n))
        mpn_neg(T[2], T[2], 2*n);
    // T3 = T2^2
    mpn_sqr(T[3], T[2], 2*n);
    // Rz = P_Qx*T3
    mpn_mul(T[5], T[3], 4*n, P_Qx, n);
    mpn_tdiv_qr(T[6], Rz, 0, T[5], 5*n, curve_p, n);
}

void ecc_xdbl(
    mp_limb_t* Rx, mp_limb_t* Rz,   // Rx, Rz must have n limbs allocated
    mp_limb_t* Px, mp_limb_t* Pz,   // Px, Pz must have n limbs allocated

    mp_limb_t* curve_a,             // curve_a must have n limbs allocated
    mp_limb_t* curve_b,             // curve_b must have n limbs allocated
    mp_limb_t* curve_p,             // curve_p must have n limbs allocated
    mp_size_t n,                    // number of limbs in curve->p
    
    ecc_xtemp T                     // temporary variables, allocated with ecc_init_xtemp(T, n).
)
{
    // T0 = Px^2
    mpn_sqr(T[0], Px, n);
    // T1 = Pz^2
    mpn_sqr(T[1], Pz, n);
    // T5 = a*T1
    mpn_mul(T[5], T[1], 2*n, curve_a, n);
    // T3 = |T0-T5|
    if (mpn_sub(T[3], T[5], 3*n, T[0], 2*n))
        mpn_neg(T[3], T[3], 3*n);
    // T4 = T3^2
    mpn_sqr(T[4], T[3], 3*n);
    // T3 = b*T1
    mpn_mul(T[3], T[1], 2*n, curve_b, n);
    // T2 = Px*Pz
    mpn_mul_n(T[2], Px, Pz, n);
    // T6 = T3*T2
    mpn_mul(T[6], T[3], 3*n, T[2], 2*n);
    // T6 = 8*T6
    T[6][5*n] = mpn_lshift(T[6], T[6], 5*n, 3);
    // T4 = |T4-T6|
    if (mpn_sub(T[4], T[4], 6*n, T[6], 5*n+1) == 0) {
        // T4 >= T6: Rx = T4 mod p
        mpn_tdiv_qr(T[6], Rx, 0, T[4], 6*n, curve_p, n);
    } else {
        // T4 < T6: Rx = p - T4 mod p
        mpn_neg(T[4], T[4], 6*n);
        mpn_tdiv_qr(T[6], Rx, 0, T[4], 6*n, curve_p, n);
        mpn_sub_n(Rx, curve_p, Rx, n);
    }
    // T6 = T0+T5
    T[6][3*n] = mpn_add(T[6], T[5], 3*n, T[0], 2*n);
    // T0 = T2*T6
    mpn_mul(T[0], T[6], 3*n+1, T[2], 2*n);
    // T5 = T3*T1
    mpn_mul(T[5], T[3], 3*n, T[1], 2*n);
    // T0 = T0+T5
    T[0][5*n+1] = mpn_add(T[0], T[0], 5*n+1, T[5], 5*n);
    // Rz = 4*T0
    mpn_lshift(T[0], T[0], 5*n+2, 2);                       // Because highest limb == 0 or 1, so shift 2 probably doesn't overflow...
    mpn_tdiv_qr(T[4], Rz, 0, T[0], 5*n+2, curve_p, n);
}

void ecc_xmul(
    mp_limb_t* Rx, mp_limb_t* Rz,       // Rx, Rz must have n limbs allocated
    mp_limb_t* Px, mp_limb_t* Pz,       // Px, Pz must have n limbs allocated

    mp_limb_t* k,                       // k must have n limbs allocated

    mp_limb_t* curve_a,                 // curve_a must have n limbs allocated
    mp_limb_t* curve_b,                 // curve_b must have n limbs allocated
    mp_limb_t* curve_p,                 // curve_p must have n limbs allocated
    mp_size_t n,                        // number of limbs in curve->p
    
    ecc_xtemp T                         // temporary variables, allocated with ecc_init_xtemp(T, n).
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
    mpn_copyi(Rz, Pz, n);

    // R1 = P*2
    mp_limb_t* Tx = T[7];
    mp_limb_t* Tz = T[8];
    ecc_xdbl(
        Tx, Tz,
        Px, Pz,

        curve_a,
        curve_b,
        curve_p,
        n,

        T
    );

    /* ---------- Doing ladder ---------- */
    while (i>=0) {
        k_limb = k[i];
        while (j>=0) {
            if ((k_limb >> j) & 1) {
                // R0 <- R1 + R0
                ecc_xadd(
                    Rx, Rz,
                    Rx, Rz,
                    Tx, Tz,
                    Px, Pz,

                    curve_a,
                    curve_b,
                    curve_p,
                    n,

                    T
                );

                // R1 <- 2*R1
                ecc_xdbl(
                    Tx, Tz,
                    Tx, Tz,

                    curve_a,
                    curve_b,
                    curve_p,
                    n,

                    T
                );
            }
            else {
                // R1 <- R1 + R0
                ecc_xadd(
                    Tx, Tz,
                    Tx, Tz,
                    Rx, Rz,
                    Px, Pz,

                    curve_a,
                    curve_b,
                    curve_p,
                    n,

                    T
                );

                // R0 <- 2*R0
                ecc_xdbl(
                    Rx, Rz,
                    Rx, Rz,

                    curve_a,
                    curve_b,
                    curve_p,
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

int ecc_xz_to_X(
    mp_limb_t* PX,                      // PX must have n limbs allocated
    mp_limb_t* Px, mp_limb_t* Pz,       // Px, Pz must have n limbs allocated

    mp_limb_t* curve_p,                 // curve_p must have n limbs allocated
    mp_size_t n,
    
    ecc_xtemp T                         // temporary variables, allocated with ecc_init_xtemp(T, n).
)
{
    mpn_copyd(T[0], Pz, n);
    mpn_copyd(T[1], curve_p, n);

    mp_size_t sn;
    mpn_zero(T[6], n);
    mpn_gcdext(T[5], T[6], &sn, T[0], n, T[1], n);

    // curve_p divides Pz: return no.
    if (!sn)
        return 0;

    // negative T[6]: T[6] = p - T[6]
    if (sn < 0)
        mpn_sub_n(T[6], curve_p, T[6], n);

    // PX = Px / Pz mod p
    mpn_mul_n(T[2], Px, T[6], n);
    mpn_tdiv_qr(T[7], PX, 0, T[2], 2*n, curve_p, n);

    // todo: i hope to remove this shit
    if (abs(sn) > n) { 
        printf("[error] wtf, we have to take care of this shit, at ecc_xz_to_X() where inverted value > p????\n");
        exit(-1);
    }

    return 1;
}