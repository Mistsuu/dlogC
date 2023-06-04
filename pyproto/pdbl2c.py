import sys

filename = sys.argv[1]
lines = open(filename, "r").readlines()
lines = [line.strip() for line in lines]

py2cmap = {
    "X1": "Px",
    "Y1": "Py",
    "Z1": "Pz",
    "X2": "Qx",
    "Y2": "Qy",
    "Z2": "Qz",
    "X3": "Rx",
    "Y3": "Ry",
    "Z3": "Rz",
    "a": "curve_aR",
    "b": "curve_bR",
    "b3": "curve_b3R",

    "XX": "T[0]",
    "ZZ": "T[1]",
    "t0": "T[2]",
    "w": "T[3]",
    "s": "T[4]",
    "ss": "T[5]",
    "R": "T[6]",
    "RR": "T[7]",
    "t3": "T[8]",
    "B": "T[9]",
    "h": "T[10]",
}

# mul_overlap_replace = 'T[6]'
# assert mul_overlap_replace not in list(py2cmap.values())

mul_temp_replace = 'T[11]'
assert mul_temp_replace not in list(py2cmap.values())

def map_to_Cvar(token):
    if token in py2cmap:
        return py2cmap[token]
    print(f'[warning] variable "{token}" does not have a corresponding map!')
    return token

result = []
for line in lines:
    result.append(f'// {line}\n')
    var_left, expr_right = line.split(' = ')
    op3 = map_to_Cvar(var_left)

    if len(tokens := expr_right.split('+')) == 2:
        op1 = map_to_Cvar(tokens[0])
        op2 = map_to_Cvar(tokens[-1])

        if op1 != op2:
            result.append(
                f'mpn_montgomery_addmod_n({op3}, {op1}, {op2}, curve_p, n);\n'
            )
        else:
            result.append(
                f'mpn_montgomery_lshift1mod_n({op3}, {op1}, curve_p, n);\n'
            )
        continue

    if len(tokens := expr_right.split('-')) == 2:
        op1 = map_to_Cvar(tokens[0])
        op2 = map_to_Cvar(tokens[-1])
        result.append(
            f'mpn_montgomery_submod_n({op3}, {op1}, {op2}, curve_p, n);\n'
        )
        continue

    if len(tokens := expr_right.split('*')) == 2:
        op1 = map_to_Cvar(tokens[0])
        op2 = map_to_Cvar(tokens[-1])
        # void mpn_montgomery_mulmod_n(
        #     mp_limb_t* rp, 
        #     const mp_limb_t* s1p, const mp_limb_t* s2p, 
        #     const mp_limb_t* dp, const mp_limb_t* Dp, 
        #     mp_size_t n, 
        #     mp_limb_t* tp
        # );

        # if multiply does not allow overlap
        # use this.
        # if op3 == op1 or op3 == op2:
        #     if var_left in py2cmap:
        #         py2cmap[var_left], mul_overlap_replace = mul_overlap_replace, op3
        #         op3 = py2cmap[var_left]
        #     else:
        #         result.extend([
        #             f'mpn_montgomery_mulmod_n({mul_overlap_replace}, {op1}, {op2}, curve_p, curve_P, n, {mul_temp_replace});\n',
        #             f'mpn_copyd({op3}, {mul_overlap_replace}, n);\n',
        #         ])
        #         continue

        result.append(
            f'mpn_montgomery_mulmod_n({op3}, {op1}, {op2}, curve_p, curve_P, n, {mul_temp_replace});\n'
        )
        continue

    if len(tokens := expr_right.split('*')) == 3:
        op1 = map_to_Cvar(tokens[0])
        op2 = map_to_Cvar(tokens[-1])
        assert op2 == '2' and tokens[1] == ''

        # void mpn_montgomery_sqrmod_n(
        #     mp_limb_t* rp, 
        #     const mp_limb_t* s1p, 
        #     const mp_limb_t* dp, const mp_limb_t* Dp,
        #     mp_size_t n, 
        #     mp_limb_t* tp
        # );

        result.append(
            f'mpn_montgomery_sqrmod_n({op3}, {op1}, curve_p, curve_P, n, {mul_temp_replace});\n'
        )
        continue
    
open(filename.split('.')[0] + '2c.c', "w").writelines(result)