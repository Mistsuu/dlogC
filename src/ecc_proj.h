#ifndef ECC_PROJ_H
#define ECC_PROJ_H

#include <gmp.h>

typedef mp_limb_t* ecc_ptemp[15];

void ecc_init_ptemp(ecc_ptemp T, mp_size_t n);
void ecc_free_ptemp(ecc_ptemp T);
void ecc_padd(
    mp_limb_t* Rx,  mp_limb_t* Ry,  mp_limb_t* Rz,       // Rx, Ry, Rz must have n limbs allocated
    mp_limb_t* Px,  mp_limb_t* Py,  mp_limb_t* Pz,       // Px, Py, Pz must have n limbs allocated
    mp_limb_t* QxR, mp_limb_t* QyR,                      // QxR, QyR   must have n limbs allocated

    mp_limb_t* curve_p,                                  // curve_p must have n limbs allocated
    mp_limb_t* curve_P,                                  // curve_P must have n limbs allocated
    mp_size_t n,                                         // number of limbs in curve->p
    
    ecc_ptemp T                                          // temporary variables, allocated with ecc_init_ptemp(T, n).
);

void ecc_pdbl(
    mp_limb_t* Rx, mp_limb_t* Ry, mp_limb_t* Rz,         // Rx, Ry, Rz must have n limbs allocated
    mp_limb_t* Px, mp_limb_t* Py, mp_limb_t* Pz,         // Px, Py, Pz must have n limbs allocated

    mp_limb_t* curve_aR,                                 // curve_aR must have n limbs allocated
    mp_limb_t* curve_p,                                  // curve_p must have n limbs allocated
    mp_limb_t* curve_P,                                  // curve_P must have n limbs allocated
    mp_size_t n,                                         // number of limbs in curve->p
    
    ecc_ptemp T                                          // temporary variables, allocated with ecc_init_ptemp(T, n).
);

int ecc_peqx(
    mp_limb_t* Px, mp_limb_t* Pz,       // Px, Pz must have n limbs allocated
    mp_limb_t* Qx, mp_limb_t* Qz,       // Qx, Qz must have n limbs allocated

    mp_limb_t* curve_p,                 // curve_p must have n limbs allocated
    mp_limb_t* curve_P,                 // curve_P must have n limbs allocated
    mp_size_t n,                        // number of limbs in curve->p
    
    ecc_ptemp T                         // temporary variables, allocated with ecc_init_ptemp(T, n).
);

int ecc_peq(
    mp_limb_t* Px, mp_limb_t* Py, mp_limb_t* Pz,  // Px, Py, Pz must have n limbs allocated
    mp_limb_t* QxR, mp_limb_t* QyR,               // QxR, QyR   must have n limbs allocated

    mp_limb_t* curve_p,                           // curve_p must have n limbs allocated
    mp_limb_t* curve_P,                           // curve_P must have n limbs allocated
    mp_size_t n,                                  // number of limbs in curve->p
    
    ecc_ptemp T                                   // temporary variables, allocated with ecc_init_ptemp(T, n).
);

void ecc_pxz_to_X(
    mp_limb_t* PX,                      // PX must have n limbs allocated
    mp_limb_t* Px, mp_limb_t* Pz,       // Px, Pz must have n limbs allocated

    mp_limb_t* curve_p,                 // curve_p must have n limbs allocated
    mp_size_t n,
    
    ecc_ptemp T                         // temporary variables, allocated with ecc_init_ptemp(T, n).
);

#endif