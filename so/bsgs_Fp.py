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

    # Check
    if k != None:
        assert pow(G, k, p) == kG % p, ValueError(f"Error here: G={G}, kG={kG}, p={p}, k={k}")
    return k
    
if __name__ == '__main__':
    p = 17964272748708788969244458806302313036635817955311740380168420334513581834233812715299218206666201505235752340555208769745447948261331590897346586194446239409383713825475563106313546766458796807501540631034204460080562480950404777411064046623626989383367603590049699912119299791858052353966829770598609512824833
    g = 9412947043258226763029447485625324661003717591200017306608138518026370832572310787843921940883700576459709515981482491065879503620142335304238377114949351468734562973562717491957611354816470822634673015638074233271308652300680823378765215437107529135162395647546397351066519737618863869091184821118645289741682
    h = 13926804676619179835395807624042712463380567581429039907413012691109926537391586248304347130233644221004746524841187952910173027927262059806754094757416806030828583346262476728771951587715305385172726807250715772624366259096660028719223127482952786188747294426332305463388900635821567015767961890254565038938347
    g_ord = 2521604123

    # Pre-calculated result
    x = 167594669
    assert pow(g, x, p) == h % p
    
    # Using algorithm to find out :)
    recoveredX = discrete_log_babystep_giantstep_Fp(
                    g, 
                    h, 
                    p, 
                    g_ord,
                    ncores=1
                )
    print(f'{recoveredX = }')
    assert pow(g, recoveredX, p) == h % p