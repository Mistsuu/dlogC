#   
#   This formula is not made by me! It is based on this website that you can
#   look for similar ones:
#       https://hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#doubling-dbl-2015-rcb
#       https://hyperelliptic.org/EFD/g1p/auto-shortw-projective.html
#    

def xDBL(X1, Y1, Z1, a, b3):
    t0 = X1**2
    t1 = Y1**2
    t2 = Z1**2
    t3 = X1*Y1
    t3 = t3+t3
    Z3 = X1*Z1
    Z3 = Z3+Z3
    X3 = a*Z3
    Y3 = b3*t2
    Y3 = X3+Y3
    X3 = t1-Y3
    Y3 = t1+Y3
    Y3 = X3*Y3
    X3 = t3*X3
    Z3 = b3*Z3
    t2 = a*t2
    t3 = t0-t2
    t3 = a*t3
    t3 = t3+Z3
    Z3 = t0+t0
    t0 = Z3+t0
    t0 = t0+t2
    t0 = t0*t3
    Y3 = Y3+t0
    t2 = Y1*Z1
    t2 = t2+t2
    t0 = t2*t3
    X3 = X3-t0
    Z3 = t2*t1
    Z3 = Z3+Z3
    Z3 = Z3+Z3