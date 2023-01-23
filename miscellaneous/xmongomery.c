#include <gmp.h>
#include <string.h>

/*
    ecc_xadd: 
        ? Add 2 points with coordinate (X:Z) form on elliptic curve.
        ! No overlap should occur between PQx, PQz and Rx, Rz or Px, Pz.
*/
void ecc_xadd(
    mp_limb_t* Rx, mp_limb_t* Rz,   // Rx, Rz must have n limbs allocated.
    mp_limb_t* Px, mp_limb_t* Pz,   // Px, Pz must have n limbs allocated
    mp_limb_t* Qx, mp_limb_t* Qz,   // Qx, Qz must have n limbs allocated
    mp_limb_t* PQx, mp_limb_t* PQz, // PQx, PQz must have n limbs allocated (P - Q)

    mp_limb_t* curve_p,             // curve_p must have n limbs allocated
    mp_size_t n,                    // number of limbs in curve->p

    // Temporary variables.
    // V[0] must have n+1 limbs allocated
    // V[1] must have 2n+2 limbs allocated
    // V[2] must have 2n+2 limbs allocated
    // V[3] must have 5n+6 limbs allocated
    // V[4] must have 5n+6 limbs allocated
    mp_limb_t* V[5] 
)
{
    /* V0 = Px + Pz */
    V[0][n] = mpn_add_n(V[0], Px, Pz, n);               // todo: check if this carry is correct.
    
    /* V1 = Qx - Qz */
    V[1][n] = mpn_sub_n(V[1], Qx, Qz, n);                     
    if (V[1][n]) 
        mpn_add_n(V[1], V[1], curve_p, n);
    V[1][n] ^= V[1][n];

    /* V1 = V1 * V0 */
    mpn_mul_n(V[2], V[1], V[0], n+1);                         
    mpn_copyi(V[1], V[2], 2*n+2);                       // todo: check if this copy is correct.
    
    /* V0 = Px - Pz */                                               
    V[0][n] = mpn_sub_n(V[0], Px, Pz, n);
    if (V[0][n])
        mpn_add_n(V[0], V[0], curve_p, n);
    V[0][n] ^= V[0][n];

    /* V2 = Qx + Qz */
    V[2][n] = mpn_add_n(V[2], Qx, Qz, n);               // todo: check if this carry is correct.

    /* V2 = V2 * V0 */
    mpn_mul_n(V[3], V[2], V[0], n+1);                         
    mpn_copyi(V[2], V[3], 2*n+2);                       // todo: check if this copy is correct.

    /* V3 = V1 + V2 */         
    V[3][2*n+2] = mpn_add_n(V[3], V[1], V[2], 2*n+2);   // todo: check if this carry is correct.

    /* V4 = V3 * V3 */
    mpn_sqr(V[4], V[3], 2*n+3);

    /* V3 = PQz * V4 */
    mpn_mul(V[3], PQz, n, V[4], 4*n+6);

    /* Rx = V3 mod p */
    mpn_tdiv_qr(V[4], Rx, 0, V[3], 5*n+6, curve_p, n);

    /* V4 = |V1 - V2| */
    if (mpn_cmp(V[1], V[2], 2*n+2) < 0) mpn_sub_n(V[4], V[2], V[1], 2*n+2);
    else                                mpn_sub_n(V[4], V[1], V[2], 2*n+2);

    /* V3 = V4 * V4 */
    mpn_sqr(V[3], V[4], 2*n+2);

    /* V4 = PQx * V3 */
    mpn_mul(V[4], PQx, n, V[3], 4*n+4);

    /* Rz = V4 mod p */
    mpn_tdiv_qr(V[3], Rz, 0, V[4], 5*n+4, curve_p, n);
}


/*
    ecc_xadd: 
        ? Double a point coordinate (X:Z) form on elliptic curve.
*/
void ecc_xdbl(
    mp_limb_t* Rx, mp_limb_t* Rz,   // Rx, Rz must have n limbs allocated
    mp_limb_t* Px, mp_limb_t* Pz,   // Px, Pz must have n limbs allocated

    mp_limb_t* curve_p,             // curve_p must have n limbs allocated
    mp_limb_t* curve_A,             // curve_A must have n limbs allocated (curve_A = (curve->a + 2) / 4) 
    mp_size_t n,                    // number of limbs in curve->p
    
    // Temporary variables.
    // V[0] must have n+1 limbs allocated
    // V[1] must have 2n+2 limbs allocated
    // V[2] must have 2n+2 limbs allocated
    // V[3] must have 5n+6 limbs allocated
    // V[4] must have 5n+6 limbs allocated
    mp_limb_t* V[5]
)
{
    /* V0 = Px + Pz; */
    V[0][n] = mpn_add_n(V[0], Px, Pz, n);

    /* V1 = V0 * V0; */
    mpn_sqr(V[1], V[0], n+1);

    /* V0 = Px - Pz; */
    V[0][n] = mpn_sub_n(V[0], Px, Pz, n);
    if (V[0][n])
        mpn_add_n(V[0], V[0], curve_p, n);
    V[0][n] ^= V[0][n];

    /* V2 = V0 * V0; */
    mpn_sqr(V[2], V[0], n+1);

    /* V3 = V1 * V2 */
    mpn_mul_n(V[3], V[1], V[2], 2*n+2);

    /* Rx = V3 mod p */
    mpn_tdiv_qr(V[4], Rx, 0, V[3], 4*n+4, curve_p, n);

    /* V1 = (V1 - V2) mod p */
    if (mpn_cmp(V[1], V[2], 2*n+2) < 0) {
        mpn_sub_n(V[1], V[2], V[1], 2*n+2);
        mpn_tdiv_qr(V[3], V[1], 0, V[1], 2*n+2, curve_p, n);
        mpn_sub_n(V[1], curve_p, V[1], n);
        mpn_zero(&V[1][n], n+2);
    } 
    else {
        mpn_sub_n(V[1], V[1], V[2], 2*n+2);
    }

    /* V3 = curve_A * V1; */
    mpn_mul(V[3], V[1], 2*n+2, curve_A, n);

    /* V3 = V3 + V2; */
    V[3][3*n+2] = mpn_add(V[3], V[3], 3*n+2, V[2], 2*n+2);

    /* V4 = V1 * V3 */
    mpn_mul(V[4], V[1], 2*n+2, V[3], 3*n+3);
    
    /* Rz = V4 mod p */
    mpn_tdiv_qr(V[3], Rz, 0, V[4], 5*n+5, curve_p, n);
}

void ecc_xmul(
    mp_limb_t* Rx, mp_limb_t* Rz,   // Rx, Rz must have n limbs allocated
    mp_limb_t* Px, mp_limb_t* Pz,   // Px, Pz must have n limbs allocated

    mp_limb_t* k,                   // multiplier k, must have n limbs allocated

    mp_limb_t* curve_p,             // curve_p must have n limbs allocated
    mp_limb_t* curve_A,             // curve_A must have n limbs allocated (curve_A = (curve->a + 2) / 4) 
    mp_size_t n,                    // number of limbs in curve->p

    // Temporary variables.
    // V[0] must have n+1 limbs allocated
    // V[1] must have 2n+2 limbs allocated
    // V[2] must have 2n+2 limbs allocated
    // V[3] must have 5n+6 limbs allocated
    // V[4] must have 5n+6 limbs allocated
    // V[5] - V[6] must have n limbs allocated
    mp_limb_t* V[7]
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
    mp_limb_t* Tx = V[5];
    mp_limb_t* Tz = V[6];
    ecc_xdbl(
        Tx, Tz,
        Px, Pz,

        curve_p,
        curve_A,
        n,

        V
    );

    /* ---------- Doing Mongomery's ladder ---------- */
    while (i>=0) {
        k_limb = k[i];
        while (j>=0) {
            printf("R0: GF(p)("); mpn_printf(Rx, n); printf(")/"); mpn_printf(Rz, n); printf("\n");
            printf("R1: GF(p)("); mpn_printf(Tx, n); printf(")/"); mpn_printf(Tz, n); printf("\n");
            printf("bit: %d\n", (k_limb >> j) & 1);

            if ((k_limb >> j) & 1) {
                // R0 <- R1 + R0
                ecc_xadd(
                    Rx, Rz,
                    Rx, Rz,
                    Tx, Tz,
                    Px, Pz,

                    curve_p,
                    n,

                    V
                );

                // R1 <- 2*R1
                ecc_xdbl(
                    Tx, Tz,
                    Tx, Tz,

                    curve_p,
                    curve_A,
                    n,

                    V
                );
            }
            else {
                // R1 <- R1 + R0
                ecc_xadd(
                    Tx, Tz,
                    Tx, Tz,
                    Rx, Rz,
                    Px, Pz,

                    curve_p,
                    n,

                    V
                );

                // R0 <- 2*R0
                ecc_xdbl(
                    Rx, Rz,
                    Rx, Rz,

                    curve_p,
                    curve_A,
                    n,

                    V
                );
            }
            j--;

            printf("R0: GF(p)("); mpn_printf(Rx, n); printf(")/"); mpn_printf(Rz, n); printf("\n");
            printf("R1: GF(p)("); mpn_printf(Tx, n); printf(")/"); mpn_printf(Tz, n); printf("\n");
            printf("---------------------------------------\n");
        }
        i--;
        j = mp_bits_per_limb-1;
    }

    // R <- R0 (already done)
}