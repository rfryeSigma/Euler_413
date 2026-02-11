"""
Run a brute force search over small rationals for solutions using Pieza's pairs method.
https://math.stackexchange.com/questions/4857229/on-why-solutions-to-x4y4z4-1-come-in-pairs

"""
from solutions import known
from datetime import datetime
from itertools import combinations
from math import isqrt
from timeit import repeat, time

def abcd_to_xyz(abcd: tuple) -> tuple:
    """Convert (a, b, c, d) representing a^4 + b^4 + c^4 = d^4 
    to rational form (x, y, z) representing x^4 + y^4 + z^4 = 1
    """
    A, B, C, D = map(QQ, abcd)
    x, y, z = A/D, B/D, C/D
    return x, y, z

def xyz_to_u(x :Integer, y :Integer, z :Integer) -> Rational:
    """Convert (x, y, z) representing x^4 + y^4 + z^4 = 1
    to single parameter u.
    """
    x_y = x - y
    u = (x_y*x_y - z*z - 1) / (x*x - x*y + y*y + x_y)
    return u

def abcd_to_u_set(abcd: tuple) -> set:
    """Convert (a, b, c, d) to (x, y, z) and then consider
    2 permutations and 8 combinations of signs on (x, y, z).
    Since the u, v, w formulas rotate (x, y, z), only swap y, z.
    to generate the set of all EC parameters u.
    There will be 4 duplicates, so return unique 12.
    """
    u_set = set()
    xx, yy, zz = abcd_to_xyz(abcd)
    for x, y, z in ((xx, yy, zz), (zz, xx, yy), (zz, yy, xx),):
        for sx in (x, -x):
            for sy in (y, -y):
                for sz in (z, -z):
                    u_set.add(xyz_to_u(sx, sy, sz))
    return u_set

"""
abcd_to_u_set is faster than abcd_to_u_set_uvw
r = repeat('abcd_to_u_set_uvw(abcd)', number=100, repeat=100, globals=globals()); np.median(r); np.std(r)
np.float64(0.018145770533010364)
np.float64(0.001256843665000622)
r = repeat('abcd_to_u_set(abcd)', number=100, repeat=100, globals=globals()); np.median(r); np.std(r)
np.float64(0.009961916599422693)
np.float64(0.001336732926780821)
"""

def known_to_bounds(known_limit: int=1000, 
                    m_limit: int=int(1e6), n_limit: int=int(1e7), ) -> tuple:
    """Search known solutions to estimate bound on u, v.
    
    Each u consists of a postive m=4k and a signed odd with with abs n.
    For each u in the set of 12 u for a solution,
    select the two u with smallest m and the 2 with smallest n.
    They will probably overlap.
    Keep the two n that are paired with the smallest m, 
    and the two m that are paired with the smallest n.
    The (m, n) bounds for u are the maxima of the first values,
    and the (m, n) bounds for v are the maxima of the second values.
    """
    um = [] # m paired with smallest n
    un = [] # n paired with smallest m
    vm = [] # m paired with 2nd smallest n
    vn = [] # n paired with 2nd smallest m
    for inx_known, val in enumerate(known.values(), start=1):
        abcd = val['abcd']
        u_set = abcd_to_u_set(abcd)
        mn = [0] * 12 # (m, n) for each u to be sorted
        nm = [0] * 12 # (n, m) for each u to be sorted
        for inx_u, u in enumerate(u_set):
            u = abs(u)
            m, n = numerator(u), denominator(u)
            if m % 2 != 0: m, n = n, m
            mn[inx_u] = m, n
            nm[inx_u] = n, m
        mn.sort()
        nm.sort()
        if nm[0][1] < m_limit:
            um.append(nm[0][1])
            if nm[1][1] < m_limit:
                vm.append(nm[1][1])
            else:
                print(f'Skipping vm {nm[1][1]} in #{inx_known}')
        else:
            print(f'Skipping um {nm[0][1]} and vm in #{inx_known}')

        if mn[0][1] < n_limit:
            un.append(mn[0][1])
            if mn[1][1] < n_limit:
                vn.append(mn[1][1])
            else:  
                print(f'Skipping vn {mn[1][1]} in #{inx_known}')
        else:  
            print(f'Skipping un {mn[0][1]} and vn in #{inx_known}')
        if len(um) >= known_limit:
            break
    return (f'{max(um):_}, {max(un):_}'), (f'{max(vm):_}, {max(vn):_}')
"""
known_to_bounds()
Skipping vm 19835764 in #10
Skipping un 11846053 and vn in #10
Skipping vm 12642040 in #54
Skipping vn 57878913 in #54
Skipping vm 12642040 in #56
Skipping vn 57878913 in #56
Skipping vm 10001951064 in #57
Skipping vn 8204718073 in #57
Skipping vm 29393447736 in #58
Skipping vn 14486729065 in #58
Skipping vm 4669000304 in #65
Skipping vn 944254963 in #65
Skipping vm 6985268 in #67
Skipping vn 92654145 in #67
Skipping vm 22869016 in #69
Skipping vm 116419537680 in #72
Skipping vn 3286280601165 in #72
Skipping vm 4669000304 in #73
Skipping vn 944254963 in #73
Skipping vm 6985268 in #74
Skipping vn 92654145 in #74
Skipping vm 323752247040 in #75
Skipping vn 634023165233 in #75
Skipping vm 1551044 in #76
Skipping vm 1856708 in #79
Skipping vm 29393447736 in #80
Skipping vn 14486729065 in #80
Skipping vm 61008600 in #84
Skipping vn 68433257 in #84
Skipping vm 4372152 in #85
Skipping vm 6037149201728 in #86
Skipping vn 3384184433553 in #86
Skipping vm 107014216 in #87
Skipping vn 210232185 in #87
Skipping vn 10490417 in #88
Skipping vm 116419537680 in #90
Skipping vn 22953456067 in #90
Skipping vm 541388136 in #92
Skipping vm 3885556 in #93
Skipping vn 23685689 in #93
Skipping vm 323752247040 in #94
Skipping vn 634023165233 in #94
('175_812, 93_017', '712_772, 4_037_701')
"""

def is_square(u) -> bool:
    """Return whether a rational is a perfect square"""
    n = numerator(u)
    if 0 >= n or isqrt(n)**2 != n: return False
    d = denominator(u)
    return isqrt(d)**2 == d

def un_square(u):
    """Return isqrt of numerator and denominator of a rational"""
    n = isqrt(numerator(u))
    d = isqrt(denominator(u))
    r = QQ(n) / QQ(d)
    assert r * r == u, f'{u}, {r}'
    return r

def v_to_D_coeffs(v):
    """
    Returns pre-calculated coefficients for the polynomial in u.
    Returns a list [c0, c1, c2, c3, c4] corresponding to powers of u.
    """
    v2 = v*v
    v3 = v2*v
    v4 = v3*v
    
    # Based on the symmetry matrix we derived:
    c0 = 4*v4 - 16*v3 + 32*v - 48
    c1 = -8*v4 + 48*v3 - 64*v2 + 32*v + 32
    c2 = -16*v3 + 48*v2 - 64*v
    c3 = -4*v4 + 8*v3 - 16*v2 + 48*v - 16
    c4 = v4 - 4*v3 - 8*v + 4
    
    return [c0, c1, c2, c3, c4]

def u_to_D(u, coeffs):
    """
    Evaluates the polynomial c4*u^4 + c3*u^3 + c2*u^2 + c1*u + c0
    using Horner's Method for maximum inner-loop speed.
    """
    c0, c1, c2, c3, c4 = coeffs
    D2 = c0 + u*(c1 + u*(c2 + u*(c3 + u*c4)))
    if is_square(D2): return un_square(D2)
    else: return None

def vu_to_P0(v, u):
    """Calculate the EC parameters Pk from v, u
    """
    u2 = u * u
    v2 = v * v
    P0 = (2 + u2) * (2 + v2) * (12 - 8*u + 2*u2 - 8*v + 2*v2 + u2*v2)
    return P0

def vu_to_Pk(v, u):
    """Calculate the EC parameters Pk from v, u
    """
    u2 = u * u
    u3 = u2 * u
    v2 = v * v
    v3 = v2 * v
    P0 =(2 + u2) * (2 + v2) * (12 - 8*u + 2*u2 - 8*v + 2*v2 + u2*v2)
    P1 = (-4 + 4*u + 2*v - u2*v) * (8*u -4*u2 + 4*v -8*u*v + 2*u2*v -4*v2 + 2*v3 + u2*v3)
    P2 = 2 * (-4 + 4*u + 2*v -u2*v) * (4 -2*u -4*v +u*v2) * (-u + v)
    P3 = (4 - 2*u - 4*v + u*v2) * (4*u - 4*u2 + 2*u3 + 8*v - 8*u*v - 4*v2 + 2*u*v2 + u3*v2)
    return (P0, P1, P2, P3)

def uvD_to_xyz(u, v, D):
    """Calculate (x, y, z) from u, v
    """
    P0, P1, P2, P3 = vu_to_Pk(v, u)
    u2 = u * u
    v2 = v * v
    denom = P0 + D*D
    x1 = (u2 - 2*u + 2) * (v2 - 2*v)
    y1 = 2*(u + v - 2) * (u + v - u*v)
    z1 = (v2 - 2*v + 2) * (u2 - 2*u)
    pair = [0]*2
    for inx, sD in enumerate((D, -D)):
        x = (P1 + x1 * sD) / denom
        y = (-(P1 + P2 + P3) + y1 * sD) / denom
        z = (P3 + z1 * sD) / denom
        pair[inx] = (x, y, z)
    return pair

def brute_search(first_vm: int, last_vm: int,
                 first_vn: int, last_vn: int,
                 last_um: int, last_un: int,
                 max_d: int=int(1e27)) -> None:
    """Search for solutions with v over given range
    and u ranging from minimal value to less than v.
    The rationals v and u are constructed from 
        m a postive multiple of 4
        and a signed odd n.
        The rational can be m/n or n/m
    There are duplicates and swap duplicates, but faster not to filter
    """
    # Look up known index from denominator.
    d_to_known_inx = {val['abcd'][-1]: inx 
            for inx, val in enumerate(known.values(), start=1)}
    found_knowns = set() # indexes of known solutions found in search

    # Search loop
    assert first_vm % 4 == 0 and last_vm % 4 == 0, "m must be multiple of 4"
    assert 4 <= first_vm <= last_vm, "first_vm must be at least 4"
    assert first_vn % 2 == 1 and last_vn % 2 == 1, "n must be odd"
    assert 1 <= first_vn <= last_vn
    assert last_um % 4 == 0 and last_un % 2 == 1
    assert 4 <= last_um <= last_vm, "u must be less than v"
    assert 1 <= last_un <= last_vn, "u must be less than v"
    new_xyz = list()
    uv_count = 0
    D_count = 0
    big_count = 0
    start = time.time()
    for vm in range(first_vm, last_vm + 1, 4):
        vu_set = set() # (v, u) pairs already searched to avoid duplicates
        for vn_raw in range(first_vn, last_vn + 1, 2):
            for vns in (vn_raw, -vn_raw):
                for v in (QQ(vm) / QQ(vns), QQ(vns) / QQ(vm)):
                    coeffs = v_to_D_coeffs(v)
                    for um in range(4, min(vm, last_um) + 1, 4):
                        for un_raw in range(1, min(vn_raw, last_un) + 1, 2):
                            for uns in (un_raw, -un_raw):
                                for u in (QQ(um) / QQ(uns), QQ(uns) / QQ(um)):
                                    if u == v:  continue
                                    uv_count += 1
                                    D = u_to_D(u, coeffs)
                                    if D is not None:
                                        D_count += 1
                                        pair = uvD_to_xyz(u, v, D)
                                        for xyz in pair:
                                            d = lcm([denominator(x) for x in xyz])
                                            print(f'({v}, {u}) -> d {d}')
                                            if d >= max_d:
                                                big_count += 1
                                                print(f'\tbig {float(d):.4e}')
                                                continue
                                            if d in d_to_known_inx:
                                                inx = d_to_known_inx[d]
                                                found_knowns.add(inx)
                                                print(f'\tknown #{inx}: {xyz}')
                                                continue
                                            print('\tnew{xyz}')
                                            new_xyz.append(xyz)

    elapsed = time.time() - start
    print(f'uv_count {uv_count}, D_count {D_count}, big_count {big_count}'
          f', knowns {len(found_knowns)} : {found_knowns}')
    print(f'new_xyz {len(new_xyz)}')
    for inx, xyz in enumerate(new_xyz, start=1):
        print(f'\t{inx}: {xyz}')
    print(f'elapsed: {elapsed:.4f}s')

"""
brute_search(4, 100, 1, 101)
elapsed: 4.2316s Just counting uv
6887996

brute_search(4, 100, 1, 101)
elapsed: 68.7418s Just checking existence of D, no D found.
(6887996, 0)

brute_search(4, 1000, 1, 47, 20, 9)
elapsed: 6.3267s Just checking existence of D, 1 D found.
(2398980, 1)

brute_search(1000, 1000, 47, 47, 20, 9)
(1000/47, -9/20) -> d 422481
	known #1: (414560/422481, 95800/422481, 217519/422481)
(1000/47, -9/20) -> d 1679142729
	known #17: (-632671960/1679142729, -1670617271/1679142729, 50237800/1679142729)
uv_count 400, D_count 1, big_count 0, knowns 2 : {1, 17}
new_xyz 0
elapsed: 0.0026s

brute_search(4, 1000, 1, 47, 20, 9)
(1000/47, -9/20) -> d 422481
	known #1: (414560/422481, 95800/422481, 217519/422481)
(1000/47, -9/20) -> d 1679142729
	known #17: (-632671960/1679142729, -1670617271/1679142729, 50237800/1679142729)
uv_count 2181396, D_count 1, big_count 0, knowns 2 : {1, 17}
new_xyz 0
elapsed: 5.9163s

datetime.now(); brute_search(4, 1520, 1, 1865, 1520, 1865); datetime.now()
datetime.datetime(2026, 2, 10, 23, 45, 24, 876310)

"""