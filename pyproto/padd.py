#
#   This formula is not made by me! It is based on this website that you can
#   look for similar ones:
#       https://hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#addition-add-2015-rcb
#       https://hyperelliptic.org/EFD/g1p/auto-shortw-projective.html
#       

def xADD(X1, X2, Y1, Y2, Z1, Z2, a, b3):
    t0 = X1*X2
    t1 = Y1*Y2
    t2 = Z1*Z2
    t3 = X1+Y1
    t4 = X2+Y2
    t3 = t3*t4
    t4 = t0+t1
    t3 = t3-t4
    t4 = X1+Z1
    t5 = X2+Z2
    t4 = t4*t5
    t5 = t0+t2
    t4 = t4-t5
    t5 = Y1+Z1
    X3 = Y2+Z2
    t5 = t5*X3
    X3 = t1+t2
    t5 = t5-X3
    Z3 = a*t4 
    X3 = b3*t2
    Z3 = X3+Z3
    X3 = t1-Z3
    Z3 = t1+Z3
    Y3 = X3*Z3
    t1 = t0+t0
    t1 = t1+t0
    t2 = a*t2
    t4 = b3*t4
    t1 = t1+t2
    t2 = t0-t2
    t2 = a*t2
    t4 = t4+t2
    t0 = t1*t4
    Y3 = Y3+t0
    t0 = t5*t4
    X3 = t3*X3
    X3 = X3-t0
    t0 = t3*t1
    Z3 = t5*Z3
    Z3 = Z3+t0