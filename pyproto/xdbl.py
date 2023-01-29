#   
#   This formula is not made by me! It is based on this website that you can
#   look for similar ones:
#       https://hyperelliptic.org/EFD/g1p/auto-code/shortw/xz/doubling/dbl-2002-it-2.op3
#       https://hyperelliptic.org/EFD/g1p/auto-shortw-xz.html
#    

T1 = Px^2           # 2N
T2 = Pz^2           # 2N
T3 = a*T2           # 3N
T4 = T1-T3          # 3N -- can take absolute
T5 = T4^2           # 6N
T6 = b*T2           # 3N
T7 = Px*Pz          # 2N
T8 = T6*T7          # 5N
T9 = 8*T8           # 5N+1
Rx = T5-T9          # N -- mod n
T10 = T1+T3         # 3N+1
T11 = T7*T10        # 5N+1
T12 = T6*T2         # 5N
T13 = T11+T12       # 5N+2
Rz = 4*T13          # 5N+3 -> N -- mod n