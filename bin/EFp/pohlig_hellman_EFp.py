import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from sage.all import *
from dlog_EFp import discrete_log_EFp, ECCPoint, ECC, ECCInf
from bsgs_EFp import discrete_log_babystep_giantstep_EFp
import random

def pohlig_hellman_EFp(
    G: int, kG: int, 
    q_e: int, E, order: int, 
    ncores: int = 4,
    alpha: int = 0,
    nranditems: int = 20
):
    """
        G: base
        kG: ans
        q_e: a prime power divides p-1
        E: Sage's curve
        order: Order of the curve.

        ncores: number of cores
        alpha: hard to explain
                (should be 3*log2(q) 
                 -- just set it to 0 
                 if unsure) 
        nranditems: number of random items.
                     (just set to 20 or more 
                      if you don't know)

        Solves for X mod q_e, where 
            b^X = a mod p
    """
    # Convert types
    b = E(G)
    a = E(kG)
    q_e = int(q_e)
    order = int(order)

    curve = ECC(
        a=int(E.a4()),
        b=int(E.a6()),
        p=int(E.base_ring().characteristic())
    )

    # Sanity checks
    assert is_prime_power(q_e) 
    assert is_prime(curve.p)
    assert order * b == E(0, 1, 0)

    # Seperate q, e
    q, e = factor(q_e)[0]
    q = int(q)
    e = int(e)

    # Solve for array x
    # which X = sum(x[i] * q^i)
    x = []

    # B is placed outside so that we can cache it :)
    B = (order//q) * b
    if B.is_zero():
        B_ = ECCInf
    else:
        B_ = ECCPoint(
            x=int(B.xy()[0]),
            y=int(B.xy()[1])
        )

    for j in range(e):
        A = (order//(q**(j+1))) * a
        if A.is_zero():
            A_ = ECCInf
        else:
            A_ = ECCPoint(
                x=int(A.xy()[0]),
                y=int(A.xy()[1])
            )

        # Using pollig-hellman / baby-step-giant-step 
        # to solve k that B^k = A (mod p)
        if q < 2**32:
            i = discrete_log_babystep_giantstep_EFp(
                curve,
                B_, A_,
                q,
                ncores
            )
        else:
            i = discrete_log_EFp(
                curve,
                B_, A_, 
                q,
                ncores, 
                alpha,
                nranditems
            )
        assert i != None, ValueError("pohlig_hellman(): not found.")
        x.append(i % q)
        
        # Update a_j
        a += (-x[-1] * q**j) * b

    X = sum(x[i] * q**i for i in range(e))
    return X

if __name__ == '__main__':
    p = 0x15501e0ee881eb3c088bdf844c64c3b6dcccbb9b072435f4732c1ac7
    a = 0x09acbbdb857a8016602753206f24086e83259589560f712ddb9b1eca
    b = 0x0b0d78078b004361e38c6d7827547930f369b368537814531fa79115
    E = EllipticCurve(GF(p), [a, b])

    order = 0x15501e0ee881eb3c088bdf844c6480b472b463767c819b730b9e6467
    factors = \
        [[2395862957, 4],
         [4084055527, 3],]
    
    G = E(
        0x0d18c96f773816c7909e2b305557f2f7b32901413f26db7165324b8a,
        0x00ff8d2e3eb7d1c4e965af2516e3bc5e2be2ec18895bca60146cab6b
    )
    assert G.order() == order
    x = random.randint(2, order-2)
    H = x*G

    print(f'{E = }')
    print(f'{G = }')
    print(f'{x = }')
    print(f'{H = }')
    
    R = []
    M = []
    for q, e in factors:
        q_e = q**e
        R.append(ans := pohlig_hellman_EFp(G, H, q_e, E, order))
        M.append(q_e)
        print(f'[+] ({p}, {e}) => {ans}')
    recoveredX = int(crt(R, M))
    print(f'{recoveredX = }')