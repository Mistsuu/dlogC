import sys

filename = sys.argv[1]
lines = open(filename, "r").readlines()
lines = [line.strip() for line in lines]

py2cmap = {
    "t0": "T[0]",
    "t1": "T[1]",
    "t2": "T[2]",
    "t3": "T[3]",
    "t4": "T[4]",
    "t5": "T[5]",
    "X1": "Px",
    "Y1": "Py",
    "Z1": "Pz",
    "X2": "Qx",
    "Y2": "Qy",
    "Z2": "Qz",
    "X3": "Rx",
    "Y3": "Ry",
    "Z3": "Rz",
    "a": "curve_a",
    "b3": "curve_b3",
}
mul_overlap_replace = 'T[6]'
mul_temp_replace = 'T[7]'
assert mul_overlap_replace not in py2cmap
assert mul_temp_replace not in py2cmap

def map_to_Cvar(token):
    if token in py2cmap:
        return py2cmap[token]
    return token

result = []
for line in lines:
    result.append(f'// {line}\n')
    var_left, expr_right = line.split(' = ')
    var_left = map_to_Cvar(var_left)

    if len(tokens := expr_right.split('+')) == 2:
        op1 = map_to_Cvar(tokens[0])
        op2 = map_to_Cvar(tokens[-1])
        result.append(
            f'mpn_montgomery_addmod_n({var_left}, {op1}, {op2}, curve_p, n);\n'
        )
        continue

    if len(tokens := expr_right.split('-')) == 2:
        op1 = map_to_Cvar(tokens[0])
        op2 = map_to_Cvar(tokens[-1])
        result.append(
            f'mpn_montgomery_submod_n({var_left}, {op1}, {op2}, curve_p, n);\n'
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

        # mpn_montgomery_mulmod_n does not allow rp == s1p or s2p
        if var_left != op1:
            result.append(
                f'mpn_montgomery_mulmod_n({var_left}, {op1}, {op2}, curve_p, curve_P, n, {mul_temp_replace});\n'
            )
        else:
            result.extend([
                f'mpn_montgomery_mulmod_n({mul_overlap_replace}, {op1}, {op2}, curve_p, curve_P, n, {mul_temp_replace});\n',
                f'mpn_copyd({var_left}, {mul_overlap_replace}, n);\n',
            ])


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
            f'mpn_montgomery_sqrmod_n({var_left}, {op1}, curve_p, curve_P, n, {mul_temp_replace});\n'
        )
        continue
    
open(filename.split('.')[0] + '2c.c', "w").writelines(result)