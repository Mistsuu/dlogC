#ifndef ECC_X_H
#define ECC_X_H

/*
    ecc_x.h:
        Define elliptic curve and its arithmetics using (X/Z) coordinate
        system.
*/
#include <gmp.h>
#include <stdlib.h>
#include <assert.h>

typedef mp_limb_t* ecc_xtemp[9];

void ecc_init_xtemp(ecc_xtemp T, mp_size_t n);
void ecc_free_xtemp(ecc_xtemp T);
void ecc_xadd(
    mp_limb_t* Rx, mp_limb_t* Rz,       // Rx, Rz must have n limbs allocated
    mp_limb_t* Px, mp_limb_t* Pz,       // Px, Pz must have n limbs allocated
    mp_limb_t* Qx, mp_limb_t* Qz,       // Qx, Qz must have n limbs allocated
    mp_limb_t* P_Qx, mp_limb_t* P_Qz,   // P_Qx, P_Qz must have n limbs allocated

    mp_limb_t* curve_a,                 // curve_a must have n limbs allocated
    mp_limb_t* curve_b,                 // curve_b must have n limbs allocated
    mp_limb_t* curve_p,                 // curve_p must have n limbs allocated
    mp_size_t n,                        // number of limbs in curve->p
    
    ecc_xtemp T                         // temporary variables.
);
void ecc_xdbl(
    mp_limb_t* Rx, mp_limb_t* Rz,   // Rx, Rz must have n limbs allocated
    mp_limb_t* Px, mp_limb_t* Pz,   // Px, Pz must have n limbs allocated

    mp_limb_t* curve_a,             // curve_a must have n limbs allocated
    mp_limb_t* curve_b,             // curve_b must have n limbs allocated
    mp_limb_t* curve_p,             // curve_p must have n limbs allocated
    mp_size_t n,                    // number of limbs in curve->p
    
    ecc_xtemp T                     // temporary variables.
);
void ecc_xmul(
    mp_limb_t* Rx, mp_limb_t* Rz,       // Rx, Rz must have n limbs allocated
    mp_limb_t* Px, mp_limb_t* Pz,       // Px, Pz must have n limbs allocated

    mp_limb_t* k,                       // k must have n limbs allocated
    
    mp_limb_t* curve_a,                 // curve_a must have n limbs allocated
    mp_limb_t* curve_b,                 // curve_b must have n limbs allocated
    mp_limb_t* curve_p,                 // curve_p must have n limbs allocated
    mp_size_t n,                        // number of limbs in curve->p
    
    ecc_xtemp T                         // temporary variables.
);

#endif