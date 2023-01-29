#
#   This formula is not made by me! It is based on this website that you can
#   look for similar ones:
#       https://hyperelliptic.org/EFD/g1p/auto-code/shortw/xz/diffadd/dadd-2002-it-3.op3
#       https://hyperelliptic.org/EFD/g1p/auto-shortw-xz.html
#       

def xADD(Px, Pz, Qx, Qz, P_Qx, P_Qz, a, b):
    T1 = Px*Qx          # 2N
    T2 = Pz*Qz          # 2N
    T3 = Px*Qz          # 2N
    T4 = Pz*Qx          # 2N
    T5 = a*T2           # 3N
    T6 = T1-T5          # 3N -- can take absolute
    T7 = T6**2          # 6N
    T8 = b*T2           # 3N
    T9 = 4*T8           # 3N+1
    T10 = T3+T4         # 2N+1
    T11 = T9*T10        # 5N+2
    T12 = T7-T11        # 6N -> N -- mod n
    Rx = P_Qz*T12       # 2N -> N -- mod n
    T13 = T3-T4         # 2N -- can take absolute
    T14 = T13**2        # 4N
    Rz = P_Qx*T14       # 5N -> N -- mod n
    return (Rx, Rz)