#   
#   Optimized code from xdbl.py -- reduce the number of variables required.
#   This formula is not made by me! It is based on this website that you can
#   look for similar ones:
#       https://hyperelliptic.org/EFD/g1p/auto-code/shortw/xz/doubling/dbl-2002-it-2.op3
#       https://hyperelliptic.org/EFD/g1p/auto-shortw-xz.html
#    

def xDBL(Px, Pz, a, b):
    T0 = Px**2          # 2N
    T1 = Pz**2          # 2N
    T5 = a*T1           # 3N
    T3 = T0-T5          # 3N -- can take absolute
    T4 = T3**2          # 6N
    T3 = b*T1           # 3N
    T2 = Px*Pz          # 2N
    T6 = T3*T2          # 5N
    T6 = 8*T6           # 5N+1
    Rx = T4-T6          # N -- mod n
    T6 = T0+T5          # 3N+1
    T0 = T2*T6          # 5N+1
    T5 = T3*T1          # 5N
    T0 = T0+T5          # 5N+2
    Rz = 4*T0           # 5N+2 -> N -- mod n
    return (Rx, Rz)