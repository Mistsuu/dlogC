import os
assert os.path.exists('./libbsgsfp.so'), \
    ValueError("You need to compile libbsgsfp.so first!")

import ctypes
libbsgsfp = ctypes.CDLL('./libbsgsfp.so')

from concurrent.futures import ProcessPoolExecutor

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

def _discrete_log_babystep_giantstep_Fp(
    G: int, kG: int, 
    p: int, upper_k: int, 
    ncores: int = 4,
    mem_limit: int = None
):
    # ============================== Sanity checks ==============================
    # Check if we can convert them into numbers?
    p = int(p)
    G = int(G) % p
    kG = int(kG) % p
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

def discrete_log_babystep_giantstep_Fp(
    G: int, kG: int, 
    p: int, upper_k: int, 
    ncores: int = 4,
    mem_limit: int = None
):
    with ProcessPoolExecutor(1) as executor:
        future = executor.submit(
                    _discrete_log_babystep_giantstep_Fp,
                    G, kG,
                    p, upper_k,
                    ncores,
                    mem_limit
                ) 
        
        return future.result()
    
if __name__ == '__main__':
    p = 32317006071311007300714876688669951960444102669715484032130345427524655138867890893197201411522913463688717960921898019494119559150490921095088152386448283120630877367300996091750197750389652106796057638384067568276792218642619756161838094338476170470581645852036305042887575891541065808607552399123930385521914333389668342420684974786564569494856176035326322058077805659331026192708460314150258592864177116725943603718461857357598351152301645904403697613233287231227125684710820209725157101726931323469678542580656697935045997268352998638215525143571875465189635916060981092608953001022308106506537956302528913389423
    g = 6105926089365190830546442344302443171907927211402399371553624069424364600958317794704621691976430060912586358762759469365450403477064599141858112932682506322699575485262882403176562256226157639166274931029457167187122430628591382618124652262598258027531381113755244608628290535616765398634271090084963303750842510271118556836687173201632007991484938122821213242223822664220175003012610505067404928235867680592563223337562759905697527412263488018223591003227557272022349488848996906639020951307104917674568523025667406667785083565411274885907958028448920971992433191149592883817115146090602878155537414327272172227013
    h = 16043821508664526231577786166314666827869569032213217280781058706332953519036585188827443833951286416350268609885445190741909040186737577720487147306708399834219618125080697212947718871693768371146951817328808879280417679410701762436062729573147242623541596150880345975714357125344492833600444984980551312803421549516014059181114813331174400896897808478108603856889960258975336041315087118069783325796807046304180948132596872928703148229309537261963848425493643411156411755017072753423243781090425601507756095555933080492152701753946534767158264936527575380896590101043073212938527188610411924127379848855626598649956
    g_ord = 42586547163467

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