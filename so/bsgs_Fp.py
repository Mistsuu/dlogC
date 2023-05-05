import os
assert os.path.exists('./libbsgsfp.so'), \
    ValueError("You need to compile libbsgsfp.so first!")

import ctypes
libbsgsfp = ctypes.CDLL('./libbsgsfp.so')

# ============================== Exposed functions ==============================
"""
int sdlog(
    char* str_p,
    char** pstr_k,
    char* str_G,
    char* str_kG,
    char* str_upper_k,
    unsigned int n_threads,
    size_t mem_limit
);
"""
sdlog = libbsgsfp.sdlog
sdlog.argtypes = [
            ctypes.c_char_p,
            ctypes.POINTER(ctypes.c_char_p),
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
sdlog_free = libbsgsfp.sdlog_free
sdlog_free.argtypes = [ ctypes.c_char_p ]
sdlog_free.restype = ctypes.c_void_p


# ============================== API Call ==============================

def discrete_log_babystep_giantstep_Fp(
    G: int, kG: int, 
    p: int, upper_k: int, 
    ncores: int = 4,
    mem_limit: int = None
):
    # ============================== Sanity checks ==============================
    # Check if we can convert them into numbers?
    G = int(G)
    kG = int(kG)
    p = int(p)
    upper_k = int(upper_k)
    ncores = int(ncores)
    if mem_limit == None:
        mem_limit = 0
    else:
        mem_limit = int(mem_limit)

    # ========================== Parse & run C functions ========================
    str_p = str(p).encode() + b'\0'
    str_k = ctypes.c_char_p()
    str_G = str(G).encode() + b'\0'
    str_kG = str(kG).encode() + b'\0'
    str_upper_k = str(upper_k).encode() + b'\0'

    sdlog(
        str_p,
        ctypes.byref(str_k),
        str_G,
        str_kG,
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
    
if __name__ == '__main__':
    p = 3229626263
    k = 756602903
    g = 5
    h = 2541821831
    g_ord = 3229626262

    print(
        discrete_log_babystep_giantstep_Fp(
            g, 
            h, 
            p, 
            g_ord
        )
    )