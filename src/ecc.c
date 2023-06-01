#include "ecc.h"
#include "num.h"
#include "const.h"
#include "calc.h"
#include "rand.h"

// ---------------------------------------- ecc ------------------------------------------

void ecc_init(ecc curve, const char* a_str, const char* b_str, const char* p_str)
{
    mpz_init(curve->a);
    mpz_init(curve->b);
    mpz_init(curve->p);

    mpz_set_str(curve->a, a_str, 10);
    mpz_set_str(curve->b, b_str, 10);
    mpz_set_str(curve->p, p_str, 10);

    // 40 checks for prime should be enough?
    assert(mpz_probab_prime_p(curve->p, BASESIZE_PRIME_CHECKER));

    // Reduce a, b mod p so they are not unecessarily big
    mpz_mod(curve->a, curve->a, curve->p);
    mpz_mod(curve->b, curve->b, curve->p);
}

void ecc_printf(ecc curve)
{
    printf("Elliptic Curve y^2 = x^3 + ");
    mpz_out_str(stdout, 10, curve->a);
    printf("*x + ");
    mpz_out_str(stdout, 10, curve->b);
    printf(" in GF(");
    mpz_out_str(stdout, 10, curve->p);
    printf(")");
}

void ecc_free(ecc curve)
{
    mpz_clear(curve->a);
    mpz_clear(curve->b);
    mpz_clear(curve->p);
}

// ---------------------------------------- eccpt ------------------------------------------

void ecc_init_pt(eccpt point)
{
    mpz_init(point->x);
    mpz_init(point->y);
    mpz_init(point->z);
}

void ecc_init_pt_str(ecc curve, eccpt point, const char* x_str, const char* y_str, const char* z_str)
{
    mpz_init(point->x);
    mpz_init(point->y);
    mpz_init(point->z);

    // X-coordinate must be presented.
    assert (x_str);
    mpz_set_str(point->x, x_str, 10);
    mpz_mod(point->x, point->x, curve->p);

    // Just set y.
    if (y_str) {
        mpz_set_str(point->y, y_str, 10);
        mpz_mod(point->y, point->y, curve->p);
    }
    else {
        mpz_set_si(point->y, UNDEFINED_COORDINATE);
    }

    // Just set z.
    if (z_str) {
        mpz_set_str(point->z, z_str, 10);
        mpz_mod(point->z, point->z, curve->p);
    }
    else {
        mpz_set_si(point->z, 1);
    }

    // Verify coordinates 
    // -- could just not implement this, 
    //    but this could be useful for debug.
    int is_point_on_curve = ecc_verify_pt(curve, point);
    assert(is_point_on_curve);
}

void ecc_init_pt_pt(eccpt dst_point, eccpt src_point)
{
    mpz_init(dst_point->x);
    mpz_init(dst_point->y);
    mpz_init(dst_point->z);
    mpz_set(dst_point->x, src_point->x);
    mpz_set(dst_point->y, src_point->y);
    mpz_set(dst_point->z, src_point->z);
}

void ecc_set_pt(eccpt dst_point, eccpt src_point)
{
    mpz_set(dst_point->x, src_point->x);
    mpz_set(dst_point->y, src_point->y);
    mpz_set(dst_point->z, src_point->z);
}

void ecc_set_pt_inf(eccpt inf_point)
{
    mpz_set_si(inf_point->x, 0);
    mpz_set_si(inf_point->y, 1);
    mpz_set_si(inf_point->z, 0);
}

void ecc_normalize_z_pt(ecc curve, eccpt point)
{
    mpz_t z_inv;
    mpz_init(z_inv);

    int z_inv_exists = mpz_invert(z_inv, point->z, curve->p);
    if (z_inv_exists) {
        mpz_set_ui(point->z, 1);
        mpz_mul(point->x, point->x, z_inv);
        mpz_mod(point->x, point->x, curve->p);
        mpz_mul(point->y, point->y, z_inv);
        mpz_mod(point->y, point->y, curve->p);
    }

    mpz_clear(z_inv);
}

int ecc_verify_pt(ecc curve, eccpt point)
{
    // infinity point
    if (mpz_cmp_si(point->z, 0) == 0)
        return 1;

    // normalize point
    if (mpz_cmp_si(point->z, 1) != 0)
        ecc_normalize_z_pt(curve, point);

    // calc y^2 = x^3 + ax + b mod p
    mpz_t point_y2;
    mpz_init(point_y2);
    mpz_mul(point_y2, point->x, point->x); // x^2
    mpz_add(point_y2, point_y2, curve->a); // x^2 + a
    mpz_mul(point_y2, point_y2, point->x); // x^3 + ax
    mpz_add(point_y2, point_y2, curve->b); // x^3 + ax + b
    mpz_mod(point_y2, point_y2, curve->p); // x^3 + ax + b mod p

    // use legendre symbol if y not here
    // recover y if we can solve for y^2.
    int result = 0;
    if (mpz_cmp_si(point->y, UNDEFINED_COORDINATE) == 0) {
        result = mpz_legendre(point_y2, curve->p) >= 0;
        if (result)
            mpz_sqrtm(point->y, point_y2, curve->p);
    }
    else {
        mpz_t point_y2_2;
        mpz_init(point_y2_2);
        mpz_mul(point_y2_2, point->y, point->y);
        mpz_mod(point_y2_2, point_y2_2, curve->p);
        result = mpz_cmp(point_y2, point_y2_2) == 0;
        mpz_clear(point_y2_2);
    }

    mpz_clear(point_y2);
    return result;
}

void ecc_random_pt(ecc curve, eccpt point)
{
    mpz_set_ui(point->z, 1);
    mpz_set_si(point->y, UNDEFINED_COORDINATE);
    do {
        mpz_dev_urandomm(point->x, curve->p);
    } while (!ecc_verify_pt(curve, point));
}

void ecc_printf_pt(eccpt point)
{
    printf("(");
    mpz_out_str(stdout, 10, point->x);
    printf(" : ");
    mpz_cmp_si(point->y, UNDEFINED_COORDINATE) == 0 ? 
        printf("?") :
        mpz_out_str(stdout, 10, point->y);
    printf(" : ");
    mpz_out_str(stdout, 10, point->z);
    printf(")");
}

void ecc_free_pt(eccpt point)
{
    mpz_clear(point->x);
    mpz_clear(point->y);
    mpz_clear(point->z);
}

// ---------------------------------------- arithemetics ------------------------------------------

/*
    ? Calculate R = P + Q.
    ! (Required every input is already init-ed)
    ! (Only call if we don't want to check P, Q on curve.)
*/
void ecc_add_noverify(ecc curve, eccpt pointR, eccpt pointP, eccpt pointQ)
{
    // If P or Q is infinity.
    if (mpz_cmp_si(pointP->z, 0) == 0) {
        ecc_set_pt(pointR, pointQ);
        return;
    }

    if (mpz_cmp_si(pointQ->z, 0) == 0) {
        ecc_set_pt(pointR, pointP);
        return;
    }

    mpz_t h, h_den;      // h and its denominator.
    mpz_t Rx, Ry;        // To prevent output to the same variable as input.
    mpz_init(h);
    mpz_init(h_den);
    mpz_init(Rx);
    mpz_init(Ry);
    
    if (mpz_cmp(pointP->x, pointQ->x) == 0) {
        // R = infinity if P == -Q.
        if (mpz_cmp(pointP->y, pointQ->y) != 0) {
            ecc_set_pt_inf(pointR);
            return;
        }

        // P == Q
        // h = (3x_P^2 + a) / 2y_P
        mpz_set(h, pointP->x);                  // x_P
        mpz_mul(h, h, h);                       // x_P^2
        mpz_mul_si(h, h, 3);                    // 3x_P^2
        mpz_add(h, h, curve->a);                // 3x_P^2 + a

        mpz_set(h_den, pointP->y);              // y_P
        mpz_mul_si(h_den, h_den, 2);            // 2y_P
        mpz_invert(h_den, h_den, curve->p);     // 1/(2y_P) mod p
    }
    else {
        // P != Q
        // h = (y_Q - y_P) / (x_Q - x_P)
        mpz_set(h, pointQ->y);                  // y_Q
        mpz_sub(h, h, pointP->y);               // y_Q - y_P

        mpz_set(h_den, pointQ->x);              // x_Q
        mpz_sub(h_den, h_den, pointP->x);       // x_Q - x_P
        mpz_invert(h_den, h_den, curve->p);     // 1/(x_Q - x_P) mod p
    }

    mpz_mul(h, h, h_den);
    mpz_mod(h, h, curve->p);

    // x_R = h^2 - x_P - x_Q
    mpz_set(Rx, h);
    mpz_mul(Rx, Rx, Rx);
    mpz_sub(Rx, Rx, pointP->x);
    mpz_sub(Rx, Rx, pointQ->x);
    mpz_mod(Rx, Rx, curve->p);

    // y_R = h(x_P - x_R) - y_P
    mpz_set(Ry, pointP->x);
    mpz_sub(Ry, Ry, Rx);
    mpz_mul(Ry, Ry, h);
    mpz_sub(Ry, Ry, pointP->y);
    mpz_mod(Ry, Ry, curve->p);

    mpz_set(pointR->x, Rx);
    mpz_set(pointR->y, Ry);
    mpz_set_si(pointR->z, 1);

    mpz_clear(h);
    mpz_clear(h_den);
    mpz_clear(Rx);
    mpz_clear(Ry);
}

/*
    ? Calculate R = P * k.
    ! (Required every input is already init-ed)
    ! (Only call if we don't want to check P on curve.)
*/
void ecc_mul_noverify(ecc curve, eccpt pointR, eccpt pointP, mpz_t k)
{
    // Set S = O
    eccpt pointS;
    ecc_init_pt(pointS);
    ecc_set_pt_inf(pointS);

    // Set T = P
    eccpt pointT; 
    ecc_init_pt(pointT); 
    ecc_set_pt(pointT, pointP);

    mpz_t abs_k;
    mpz_init(abs_k);
    mpz_abs(abs_k, k);

    mp_bitcnt_t nbits_n = mpz_size(k) * mp_bits_per_limb;
    for (mp_bitcnt_t i = 0; i < nbits_n; ++i) {
        if (mpz_tstbit(abs_k, i))
            ecc_add_noverify(curve, pointS, pointS, pointT);
        ecc_add_noverify(curve, pointT, pointT, pointT);
    }

    // In case k is negative :)
    if (mpz_cmp_si(k, 0) < 0) 
        mpz_sub(pointS->y, curve->p, pointS->y);

    mpz_clear(abs_k);
    ecc_set_pt(pointR, pointS);
    ecc_free_pt(pointT);
    ecc_free_pt(pointS);
}

/*
    ? Calculate R = P - Q.
    ! (Required every input is already init-ed)
    ! (Only call if we don't want to check P, Q on curve.)
*/
void ecc_sub_noverify(ecc curve, eccpt pointR, eccpt pointP, eccpt pointQ)
{
    eccpt point_Q;
    ecc_init_pt_pt(point_Q, pointQ);
    mpz_sub(point_Q->y, curve->p, point_Q->y);          // invert Y coordinate
    ecc_add_noverify(curve, pointR, pointP, point_Q);
    ecc_free_pt(point_Q);
}

/*
    ? Checks if P = Q.
    ! (Required every input is already init-ed)
    ! (Only call if we don't want to check P, Q on curve.)
*/
int ecc_eq_noverify(ecc curve, eccpt pointP, eccpt pointQ)
{
    if (mpz_cmp_si(pointP->z, 1) != 0 && mpz_cmp_si(pointP->z, 0) != 0)
        ecc_normalize_z_pt(curve, pointP);
    if (mpz_cmp_si(pointQ->z, 1) != 0 && mpz_cmp_si(pointQ->z, 0) != 0)
        ecc_normalize_z_pt(curve, pointQ);
    if (mpz_cmp(pointP->z, pointQ->z) != 0)
        return 0;
    return (
        mpz_cmp(pointP->x, pointQ->x) == 0 && 
        mpz_cmp(pointP->y, pointQ->y) == 0
    );
}

/*
    ? Calculate R = P + Q.
    ! (Required every input is already init-ed)
*/
void ecc_add(ecc curve, eccpt pointR, eccpt pointP, eccpt pointQ)
{
    assert(mpz_cmp_si(pointP->y, UNDEFINED_COORDINATE) != 0);
    assert(mpz_cmp_si(pointQ->y, UNDEFINED_COORDINATE) != 0);
    assert(ecc_verify_pt(curve, pointP));
    assert(ecc_verify_pt(curve, pointQ));
    return ecc_add_noverify(curve, pointR, pointP, pointQ);
}

/*
    ? Calculate R = P - Q.
    ! (Required every input is already init-ed)
*/
void ecc_sub(ecc curve, eccpt pointR, eccpt pointP, eccpt pointQ)
{
    assert(mpz_cmp_si(pointP->y, UNDEFINED_COORDINATE) != 0);
    assert(mpz_cmp_si(pointQ->y, UNDEFINED_COORDINATE) != 0);
    assert(ecc_verify_pt(curve, pointP));
    assert(ecc_verify_pt(curve, pointQ));
    return ecc_sub_noverify(curve, pointR, pointP, pointQ);
}

/*
    ? Calculate R = P * k.
    ! (Required every input is already init-ed)
*/
void ecc_mul(ecc curve, eccpt pointR, eccpt pointP, mpz_t k)
{
    assert(mpz_cmp_si(pointP->y, UNDEFINED_COORDINATE) != 0);
    assert(ecc_verify_pt(curve, pointP));
    return ecc_mul_noverify(curve, pointR, pointP, k);
}

/*
    ? Checks if P = Q.
    ! (Required every input is already init-ed)
*/
int ecc_eq(ecc curve, eccpt pointP, eccpt pointQ)
{
    assert(ecc_verify_pt(curve, pointP));
    assert(ecc_verify_pt(curve, pointQ));
    return ecc_eq_noverify(curve, pointP, pointQ);
}

// ---------------------------------------- pairing ------------------------------------------

/*
    ? Evaluate g = gPQ(R).
    ? gPQ(R) is a function whose coefficients are based on 
    ? the coordinates of P and Q
    ? that div(gPQ) = [P] + [Q] - [P+Q] - [O]

    ? This function is created by dividing 2 smaller functions:
    ?     - num where div(num) = [P] + [Q] + [-P-Q] - 3*[O]
    ?     - den where div(den) = [P+Q] + [-P-Q] + 2*[O]

    ! (Required every input is already init-ed)
    ! (Only call if we don't want to check P, Q, R on curve.)
*/
void ecc_weil_gPQ(ecc curve, mpz_t g, eccpt pointP, eccpt pointQ, eccpt pointR, mpz_t n)
{

}

/*
    ? Evaluate f = fP(R).
    ? fP is a function whose coefficients are based on 
    ? coordinates of P
    ? that div(fP) = n[P] - n[O] 
    ? with n is a number such that nP = O.

    ? The function is evaluated based on smaller functions g, built up
    ? recursively using 2 functions whose div are:
    ?         2[kP]       - [2kP]     - [O]
    ?         [2kP] + [P] - [(2k+1)P] - [O]

    ! (Required every input is already init-ed)
    ! (Only call if we don't want to check P, R on curve.)
*/
void ecc_weil_fP(ecc curve, mpz_t f, eccpt pointP, eccpt pointR, mpz_t n)
{

}

/*
    ? Calculate Weil's pairing E = e(P, Q),
    ? where P, Q has order n.
    ! (Required every input is already init-ed)
    ! (Only call if we don't want to check P, Q on curve.)
    ! (and also P.z and Q.z == 1)
*/
void ecc_weil_pairing_noverify(ecc curve, mpz_t E, eccpt pointP, eccpt pointQ, mpz_t n)
{
    // ===============================================================
    // -- Step 1: Choose a random point S that 
    // -- S != P, S != Q and S != P-Q
    // ===============================================================
    eccpt pointS;
    eccpt pointP_Q;
    ecc_init_pt(pointS);
    ecc_init_pt(pointP_Q);
    ecc_sub_noverify(curve, pointP_Q, pointP, pointQ);

    do {
        ecc_random_pt(curve, pointS);
    } while (
        ecc_eq_noverify(curve, pointS, pointP)   == 1 ||
        ecc_eq_noverify(curve, pointS, pointQ)   == 1 ||
        ecc_eq_noverify(curve, pointS, pointP_Q) == 1
    );

    // ===============================================================
    // -- Step 2: Compute outputs of 
    // -- fP(Q + S), fP(S), fQ(-S), fQ(P - S)
    // ===============================================================

    // setup -S
    eccpt point_S;
    ecc_init_pt_pt(point_S, pointS);
    mpz_sub(point_S->y, curve->p, point_S->y);
    
    // setup P-S
    eccpt pointP_S;
    ecc_init_pt(pointP_S);
    ecc_sub_noverify(curve, pointP_S, pointP, pointS);

    // setup Q+S
    eccpt pointQS;
    ecc_init_pt(pointQS);
    ecc_add_noverify(curve, pointQS, pointQ, pointS);

    mpz_t fP__QS;
    mpz_t fP__S;
    mpz_t fQ___S;
    mpz_t fQ__P_S;
    mpz_init(fP__QS);
    mpz_init(fP__S);
    mpz_init(fQ___S);
    mpz_init(fQ__P_S);
    ecc_weil_fP(curve, fP__QS,  pointP, pointQS,  n);
    ecc_weil_fP(curve, fP__S,   pointP, pointS,   n);
    ecc_weil_fP(curve, fQ___S,  pointQ, point_S,  n);
    ecc_weil_fP(curve, fQ__P_S, pointQ, pointP_S, n);

    // ===============================================================
    // -- Step 3: Compute Weil's pairing
    // ===============================================================
    mpz_t ePQ_num;
    mpz_t ePQ_den;
    mpz_t ePQ;
    mpz_init(ePQ_num);
    mpz_init(ePQ_den);
    mpz_init(ePQ);
    mpz_mul(ePQ_num, fP__QS, fQ___S);
    mpz_mul(ePQ_den, fP__S,  fQ__P_S);
    mpz_invert(ePQ, ePQ_den, curve->p);
    mpz_mul(ePQ, ePQ, ePQ_num);
    mpz_mod(ePQ, ePQ, curve->p);

    // ===============================================================
    // -- Step 4: Clean up.
    // ===============================================================
    mpz_clear(ePQ_num);
    mpz_clear(ePQ_den);
    mpz_clear(ePQ);
    mpz_clear(fP__QS);
    mpz_clear(fP__S);
    mpz_clear(fQ___S);
    mpz_clear(fQ__P_S);
    ecc_free_pt(pointQS);
    ecc_free_pt(point_S);
    ecc_free_pt(pointP_S);
    ecc_free_pt(pointS);
    ecc_free_pt(pointP_Q);
}

/*
    ? Calculate Weil's pairing E = e(P, Q),
    ? where P, Q has order n.
    ! (Required every input is already init-ed)
*/
void ecc_weil_pairing(ecc curve, mpz_t E, eccpt pointP, eccpt pointQ, mpz_t n)
{
    assert(ecc_verify_pt(curve, pointP));
    assert(ecc_verify_pt(curve, pointQ));

    // make sure
    //     P * n = O
    //     Q * n = O
    eccpt pointT;
    ecc_init_pt(pointT);

    ecc_mul_noverify(curve, pointT, pointP, n);
    assert(mpz_cmp_si(pointT->z, 0) == 0);
    ecc_mul_noverify(curve, pointT, pointQ, n);
    assert(mpz_cmp_si(pointT->z, 0) == 0);

    ecc_free_pt(pointT);
    return ecc_weil_pairing_noverify(curve, E, pointP, pointQ, n);
}