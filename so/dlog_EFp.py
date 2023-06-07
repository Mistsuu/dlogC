import os
assert os.path.exists('./libdlogefp.so'), \
    ValueError("You need to compile libdlogefp.so first!")

import ctypes
libdlogefp = ctypes.CDLL('./libdlogefp.so')

from collections import namedtuple
from concurrent.futures import ProcessPoolExecutor

# ============================== Exposed functions ==============================
"""
int sdlog(
    // Curve parameters
    char* str_curve_a,
    char* str_curve_b,
    char* str_curve_p,

    // To be overwritten
    char** pstr_k,

    char* str_Gx,
    char* str_Gy,
    char* str_Gz,
    char* str_kGx,
    char* str_kGy,
    char* str_kGz,
    char* str_n,

    // Configs
    unsigned int n_threads,
    unsigned int n_cache_items,
    unsigned int n_rand_items
);
"""
sdlog = libdlogefp.sdlog
sdlog.argtypes = [
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.c_char_p,

            ctypes.POINTER(ctypes.c_char_p),
            
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.c_char_p,

            ctypes.c_uint,
            ctypes.c_uint,
            ctypes.c_uint,
        ]
sdlog.restype = ctypes.c_int

"""
void sdlog_free(
    char* str_k
);
"""
sdlog_free = libdlogefp.sdlog_free
sdlog_free.argtypes = [ ctypes.c_char_p ]
sdlog_free.restype = ctypes.c_void_p


# ============================== API Call ==============================

ECC      = namedtuple('ECC', ['a', 'b', 'p'])
ECCPoint = namedtuple('ECCPoint', ['x', 'y'])
ECCInf   = None

def _discrete_log_EFp(
    curve: ECC,
    G: ECCPoint | ECCInf, kG: ECCPoint | ECCInf,
    n: int, 
    ncores: int = 4,
    ncacheitems: int = 4,
    nranditems: int = 20,
):
    # ============================== Sanity checks ==============================
    # Check if we can convert them into numbers?
    curve_a = int(curve.a)
    curve_b = int(curve.b)
    curve_p = int(curve.p)
    
    if G == ECCInf:
        Gx = 0
        Gy = 0
        Gz = 0
    Gx = int(G.x) % curve_p
    Gy = int(G.y) % curve_p
    Gz = 1

    if kG == ECCInf:
        kGx = 0
        kGy = 0
        kGz = 0
    kGx = int(kG.x) % curve_p
    kGy = int(kG.y) % curve_p
    kGz = 1

    n = int(n)

    ncores = int(ncores)
    ncacheitems = int(ncacheitems)
    nranditems = int(nranditems)

    # ========================== Parse & run C functions ========================
    str_curve_a = str(curve_a).encode() + b'\0'
    str_curve_b = str(curve_b).encode() + b'\0'
    str_curve_p = str(curve_p).encode() + b'\0'

    str_k = ctypes.c_char_p()
    
    str_Gx = str(Gx).encode() + b'\0'
    str_Gy = str(Gy).encode() + b'\0'
    str_Gz = str(Gz).encode() + b'\0'
    str_kGx = str(kGx).encode() + b'\0'
    str_kGy = str(kGy).encode() + b'\0'
    str_kGz = str(kGz).encode() + b'\0'
    str_n = str(n).encode() + b'\0'

    sdlog(
        str_curve_a,
        str_curve_b,
        str_curve_p,

        ctypes.byref(str_k),
        
        str_Gx,
        str_Gy,
        str_Gz,
        str_kGx,
        str_kGy,
        str_kGz,
        str_n,

        ncores,
        ncacheitems,
        nranditems
    )

    # Parse value
    k = None
    if (str_k.value != b'None'):
        k = int(str_k.value)

    # Manually free string buffer
    sdlog_free(str_k)

    return k

def discrete_log_EFp(
    curve: ECC,
    G: ECCPoint | ECCInf, kG: ECCPoint | ECCInf,
    n: int, 
    ncores: int = 4,
    ncacheitems: int = 4,
    nranditems: int = 20
):
    with ProcessPoolExecutor(1) as executor:
        future = executor.submit(
            _discrete_log_EFp,
            curve,
            G, kG,
            n,
            ncores,
            ncacheitems,
            nranditems
        ) 
        
        return future.result()
    
if __name__ == '__main__':
    curve = ECC(
        a=0x060c94d30d478908,
        b=0x038b5c123bb48b94,
        p=0x06d8fefeca5d4ca3
    )

    G = ECCPoint(
        x=0x064c29405844b615,
        y=0x039e6125d48aac3d
    )

    kG = ECCPoint(
        x=78069856179243048,
        y=373302913259898234
    )

    print(
        discrete_log_EFp(
            curve,
            G, kG,
            0x06d8fefe8066085f,

            ncores=4,
            ncacheitems=2,
            nranditems=20    
        )
    )