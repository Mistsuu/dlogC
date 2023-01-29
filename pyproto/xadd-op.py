#   
#   Optimized code from xadd.py -- reduce the number of variables required.
#   This formula is not made by me! It is based on this website that you can
#   look for similar ones:
#       https://hyperelliptic.org/EFD/g1p/auto-code/shortw/xz/diffadd/dadd-2002-it-3.op3
#       https://hyperelliptic.org/EFD/g1p/auto-shortw-xz.html
#       

T0 = Px*Qx          # 2N
T1 = Pz*Qz          # 2N
T2 = Px*Qz          # 2N
T3 = Pz*Qx          # 2N
T4 = a*T1           # 3N
T0 = T0-T4          # 3N -- can take absolute
T4 = T0^2           # 6N
T0 = b*T1           # 3N
T0 = 4*T0           # 3N+1
T1 = T2+T3          # 2N+1
T5 = T0*T1          # 5N+2
T4 = T4-T5          # 6N -> N -- mod n
Rx = P_Qz*T4        # 2N -> N -- mod n
T2 = T2-T3          # 2N -- can take absolute
T3 = T2^2           # 4N
Rz = P_Qx*T3        # 5N -> N -- mod n