import os
assert os.path.exists('./libbsgsefp.so'), \
    ValueError("You need to compile libbsgsefp.so first!")

import ctypes
libbsgsefp = ctypes.CDLL('./libbsgsefp.so')

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
    char* str_upper_k,

    // Configs
    unsigned int n_threads,
    size_t mem_limit
);
"""
sdlog = libbsgsefp.sdlog
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
            ctypes.c_size_t
        ]
sdlog.restype = ctypes.c_int

"""
void sdlog_free(
    char* str_k
);
"""
sdlog_free = libbsgsefp.sdlog_free
sdlog_free.argtypes = [ ctypes.c_char_p ]
sdlog_free.restype = ctypes.c_void_p


# ============================== API Call ==============================

ECC      = namedtuple('ECC', ['a', 'b', 'p'])
ECCPoint = namedtuple('ECCPoint', ['x', 'y'])
ECCInf   = None

def _discrete_log_babystep_giantstep_EFp(
    curve: ECC,
    G: ECCPoint | ECCInf, kG: ECCPoint | ECCInf,
    upper_k: int, 
    ncores: int = 4,
    mem_limit: int = None
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

    upper_k = int(upper_k)

    ncores = int(ncores)
    if mem_limit == None:
        mem_limit = 0
    else:
        mem_limit = int(mem_limit)

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
    str_upper_k = str(upper_k).encode() + b'\0'

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
        str_upper_k,

        ncores,
        mem_limit
    )

    # Parse value
    k = None
    if (str_k.value != b'None'):
        k = int(str_k.value)

    # Manually free string buffer
    sdlog_free(str_k)

    return k

def discrete_log_babystep_giantstep_EFp(
    curve: ECC,
    G: ECCPoint | ECCInf, kG: ECCPoint | ECCInf,
    upper_k: int, 
    ncores: int = 4,
    mem_limit: int = None
):
    with ProcessPoolExecutor(1) as executor:
        future = executor.submit(
            _discrete_log_babystep_giantstep_EFp,
            curve,
            G, kG,
            upper_k,
            ncores,
            mem_limit
        ) 
        
        return future.result()
    
if __name__ == '__main__':
    curve = ECC(
        a=448019786388741247,
        b=544225411105375163,
        p=593010448435692503
    )

    G = ECCPoint(
        x=437471552757133390,
        y=354281835126765881
    )

    kG = ECCPoint(
        x=295738136557598210,
        y=89525692852745488
    )

    print(
        discrete_log_babystep_giantstep_EFp(
            curve,
            G, kG,
            593010448361862286,
            mem_limit=4*1024*1024*1024
        )
    )