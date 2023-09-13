import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from sage.all import factor, is_prime_power, is_prime, crt
from dlog_Fp import discrete_log_Fp
from bsgs_Fp import discrete_log_babystep_giantstep_Fp

def pohlig_hellman_Fp(
    G: int, kG: int, 
    q_e: int, p: int, 
    ncores: int = 4,
    alpha: int = 0,
    nranditems: int = 20
):
    """
        G: base
        kG: ans
        p: modulus 
        q_e: a prime power divides p-1

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
    a = int(kG)
    b = int(G)
    q_e = int(q_e)
    p = int(p)

    # Sanity checks
    assert is_prime_power(q_e) 
    assert is_prime(p)
    assert (p-1) % q_e == 0
    assert pow(a, p-1, p) == 1

    # Seperate q, e
    q, e = factor(q_e)[0]
    q = int(q)
    e = int(e)

    # Solve for array x
    # which X = sum(x[i] * q^i)
    x = []

    # B is placed outside so that we can cache it :)
    B = pow(b, (p-1)//q, p)

    for j in range(e):
        A = pow(a, (p-1)//(q**(j+1)), p)

        # Using pollard-rho / babystep-giantstep
        # to solve k that B^k = A (mod p)
        if q < 2**32:
            i = discrete_log_babystep_giantstep_Fp(
                B, A,
                p, q, 
                ncores
            )
        else:
            i = discrete_log_Fp(
                    B, A, 
                    p, q, 
                    ncores, 
                    alpha,
                    nranditems
                )

        assert i != None, ValueError("pohlig_hellman(): not found.")
        x.append(i % q)
        
        # Update a_j
        a *= pow(b, -x[-1] * q**j, p)
        a %= p

    X = sum(x[i] * q**i for i in range(e))
    return X

# if __name__ == '__main__':
#     import random
#     p = 17964272748708788969244458806302313036635817955311740380168420334513581834233812715299218206666201505235752340555208769745447948261331590897346586194446239409383713825475563106313546766458796807501540631034204460080562480950404777411064046623626989383367603590049699912119299791858052353966829770598609512824833
#     factors = \
#         [[2, 47],
#          [2521604123, 4],
#          [3364561957, 4],
#          [4157103299, 3],
#          [3432964477, 3],
#          [3756067001, 4],
#          [3478664351, 2],
#          [3281130893, 2],
#          [4217020723, 2],
#          [3992745643, 4],
#          [4166436133, 3]]
    
#     g = 3
#     x = random.randint(2, p-1)
#     h = pow(g, x, p)

#     print(f'{p = }')
#     print(f'{g = }')
#     print(f'{x = }')
#     print(f'{h = }')

#     R = []
#     M = []
#     for q, e in factors:
#         q_e = q**e
#         R.append(ans := pohlig_hellman_Fp(g, h, q_e, p))
#         M.append(q_e)
#         print(f'[+] ({q}, {e}) => {ans}')
#     recoveredX = int(crt(R, M))
#     print(f'{recoveredX = }')

if __name__ == '__main__':
    import random
    p = 140601987068392810555209463610249011351761931006927193734529493874226799329006719
    factors = \
        [[2, 1],
         [3, 1],
         [97561, 1],
         [1644797261, 1],
         [37482391229, 1],
         [875386295981, 1],
         [26892940314817, 1],
         [73251320436449, 1],
         [2259283924057529, 1]]
    
    g = 3
    x = random.randint(2, p-1)
    h = pow(g, x, p)

    print(f'{p = }')
    print(f'{g = }')
    print(f'{x = }')
    print(f'{h = }')
    
    R = []
    M = []
    for q, e in factors:
        q_e = q**e
        R.append((ans := pohlig_hellman_Fp(g, h, q_e, p, 4, 0, 100)))
        M.append(q_e)
        print(f'[+] ({q}, {e}) => {ans}')
    recoveredX = int(crt(R, M))
    print(f'{recoveredX = }')