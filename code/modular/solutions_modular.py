"""
Explore whether modular curves can help find new solutions to a^4 + b^4 + c^4 = d^4.
"""
import csv
from datetime import datetime, timedelta
from math import isqrt, gcd
from solutions import known
from sage.all import help, oo, pari, sage_eval, \
    Integer, QQ, Rational, hilbert_symbol, lcm, solve, var, \
    Curve, DiagonalQuadraticForm, Jacobian, PolynomialRing
from sage.rings.polynomial.polynomial_rational_flint import Polynomial_rational_flint
from timeit import repeat, time

def xyz_to_u(x :Rational, y :Rational, z :Rational) -> Rational:
    """Convert (x, y, z) representing x^4 + y^4 + z^4 = 1
    to single parameter u.
    """
    x_y = x - y
    u = (x_y*x_y - z*z - 1) / (x*x - x*y + y*y + x_y)
    return u

def abcd_to_xyz(abcd: tuple) -> tuple:
    """Convert (a, b, c, d) representing a^4 + b^4 + c^4 = d^4 
    to rational form (x, y, z) representing x^4 + y^4 + z^4 = 1
    """
    A, B, C, D = map(QQ, abcd)
    x, y, z = A/D, B/D, C/D
    return x, y, z

def abcd_to_h_u(abcd: tuple) -> set:
    """Convert (a, b, c, d) to (x, y, z) and then consider
    3 permutations and 8 combinations of signs on (x, y, z).
    Return sorted unique (height, u) pairs.
    """
    u_set = set()
    xx, yy, zz = abcd_to_xyz(abcd)
    for x, y, z in ((xx, yy, zz), (zz, xx, yy), (zz, yy, xx),):
        for sx in (x, -x):
            for sy in (y, -y):
                for sz in (z,):
                    u = xyz_to_u(sx, sy, sz)
                    h = u.height()
                    u_set.add((h, u))
    return sorted(u_set)

def check_elkies_rules(u: Rational) -> bool:
    """Check the u satisfies Elkies rules for infinitely many solutions."""
    m, n = u.numer(), u.denom()
    if m%2: m, n = n, m
    assert m%4 == 0 and abs(n)%4 == 1, f'{u} -> {m}, {n} fail mod 4'

    # TODO Check R and S rules.
    return True
"""
scan_known(6)
-9/20 True
AssertionError: 1000/47 -> 1000, 47 fail mod 4

but -9/20 paired with 1000/47 gives piezas in pairs solutions
coeffs = u_to_D_coeffs_int(1000, 47)
D = v_to_D_int(-9, 20, coeffs, 47)
495260031/441800
pair = uvD_to_xyz(QQ(1000/47), QQ(-9/20), D)
[(-50237800/1679142729, 1670617271/1679142729, 632671960/1679142729), 
(-217519/422481, -95800/422481, -414560/422481)]

If I said Elkies rules have to work on the inverse of m/n for u or v:
56,-165, -383021,380940 and 56,-165, 2644685,570612
would work but not quickly. And they have unique solutions.
Same with some of the 136,-133 and 400,-37 pairs.

So give up on Elkies rules. Just restrict to odd and 0 mod 4.
"""

def mn_to_xyt_conics(mn: Rational) -> tuple:
    """ Apply m/n to conic equations in x, y, t
        (2m^2+n^2)y^2 = -(6m^2-8mn+3n^2)x^2 -2(2m^2-n^2)x -2mn        (1) 
        (2m^2+n^2)t^2 = 4(2m^2-n^2)x^2      +8mnx         +(n^2-2m^2) (2)
    Return ((4 coeffs for y^2 eq 1),  (4 coeffs for t^2 eq 2))
    """
    m, n = mn.numer(), mn.denom()
    m2 = m*m
    n2 = n*n
    mn = m*n
    y2_coeffs = (2 * m2 + n2, 
                 -(6 * m2 - 8 * mn + 3 *n2),
                 -2 * (2 * m2 - n2), 
                 -2 * mn)
    t2_coeffs = (2 * m2 + n2,
                 4 * (2 * m2 - n2),
                 8 * mn,
                 n2 - 2 * m2)
    return y2_coeffs, t2_coeffs

def check_quadratic_solvability(D, delta):
    """ Check the Hilbert symbol on equation X^2 - Dy^2 - delta*z^2 = 0
    at infinite prime (Real Solvability)
    and at 2 and odd primes dividing coeffs
    """
    return hilbert_symbol(D, delta, -1) == 1 == hilbert_symbol(D, delta, 2)

def check_yt(mn: Rational) -> bool:
    """ Return whether both y^2 and t^2 are solvable and have a point
    """
    y2_coeffs, t2_coeffs = mn_to_xyt_conics(mn)
    for a0, a, b, c in (t2_coeffs, y2_coeffs):
        D = 4 * a * a0
        delta = b**2 - 4 * a * c
        if not check_quadratic_solvability(D, delta):
            return False
        Q = DiagonalQuadraticForm(QQ, [1, -D, -delta])
        try: # find a rational point (X, V, Z)
            point = Q.solve()
        except Exception as e: return False
    return True

def check_known_uv_inv(file_name: str='solutions_uv.csv') -> None:
    """Check hypothesis that inverses of a uv pair
    are always y^2 and t^2 solvable and have a point.
    """
    with open(file_name, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file)
        header = next(reader)
        for row in reader:
            u_n, u_d, v_n, v_d= [int(row[i]) for i in range(4)]
            u = QQ(u_n)/u_d
            assert check_yt(u.inverse()), f'row {row} fails u_d/u_n'
            v = QQ(v_n)/v_d
            assert check_yt(v.inverse()), f'row {row} fails v_d/v_n'
            #print(row, u, v)
""" Success
"""

def u_to_quartic(u: Rational) -> Polynomial_rational_flint:
    """ Given u, build rational quartic polynomial in v.
    """
    u2 = u * u
    u3 = u2 * u
    u4 = u2 * u2
    d0 = 4 * (-6 - 2 * u + u2) * (2 - 2 * u + u2) 
    d1 = - 8 * (-2 - 4 * u + u2) * (2 - 2 * u + u2)
    d2 = - 16 * u * (4 - 3 * u + u2)
    d3 = - 4 * (4 - 12 * u + 4 * u2 - 2 * u3 + u4) 
    d4 = (4 - 8 * u - 4 * u3 + u4)
    R = PolynomialRing(QQ, 'v, D')
    v, D = R.gens()
    D2 = d0 + d1 * v + d2 * v**2 + d3 * v**3 + d4 * v**4
    D2 = D2.univariate_polynomial()
    return D2

# unfinished
def scan_known(n_known: int=8, max_h: int=100_000_000,
               max_d: int=int(1e27)) -> None:
    """For all abcd in first n_known solutions, generate 12 u.
    For each u with height less than max_h:
        Verify no obstructions on u
    ...
    """
    for val in list(known.values())[:n_known]:
        abcd = val['abcd']
        print(abcd[-1])
        for h, u in abcd_to_h_u(abcd):
            if h > max_h: continue
            #print(u, check_elkies_rules(u))
            assert check_yt(QQ(1)/u), f'No point for u {u}'
            print('\t', u)
            #TODO generate JD2 (smarter code) for each u
            # I think I should use u, but maybe 1/u?
            # Decide by running on C0, C1, C5, C7, C9, Cx
    """ More TODO
    Search for square quartic point.
    Get conductor factored, torsion points, gens
    height_pairing_matrix, new points.
    Check whether known solutions.

    map u to E, conductor, number of solutions.
    Check whether isogenous with another curve.
    """

def report_uvw(n_known: int=8):
    """ Scan known solutions for modular patterns in u,v,w.
    """
    for val in list(known.values())[:n_known]:
        t, u, v, w = val['tuvw']
        if t%5 != 0:
            t, u = u, t
            assert t%5 == 0
        a = 8 * t
        b = 8 * u
        c = w - v
        d = w + v
        assert a**4 + b**4 + c**4 == d**4
        qd = QQ(d)
        x = a / qd
        y = b / qd
        z = c / qd
        un = (x - y)**2 - z**2 - 1
        vn = (y - z)**2 - x**2 - 1
        wn = (z - x)**2 - y**2 - 1
        ud = x**2 - x*y + y**2 + (x - y)
        vd = y**2 - y*z + z**2 + (y - z)
        wd = z**2 - z*x + x**2 + (z - x)
        u = un / ud
        v = vn / vd
        w = wn / wd
        u_num, u_den = u.numer(), u.denom()
        v_num, v_den = v.numer(), v.denom()
        w_num, w_den = w.numer(), w.denom()
        print(a, b, c, d)
        print('\tu = ', u_num, '/', u_den, '=', float(u_num)/float(u_den)) 
        print('\tv = ', v_num, '/', v_den, '=', float(v_num)/float(v_den))
        print('\tw = ', w_num, '/', w_den, '=', float(w_num)/float(w_den))

def uv_to_Pk(u: Rational, v: Rational) -> tuple:
    """Calculate the parameters Pk from u, v
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

def uvD_to_xyz(u: Rational, v: Rational, D: Rational) -> tuple:
    """Calculate (x, y, z) from u, v
    """
    P0, P1, P2, P3 = uv_to_Pk(u, v)
    denom = P0 + D*D
    u2 = u * u
    v2 = v * v
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

def uv_to_w(u: Rational, v: Rational) -> Rational:
    """The u, v, w derived from x^4 + y^4 + z^4 = 1
    have relationship 2(u+v+w)-uvw-4=0
    """
    w = 2 * (u + v - 2) / (u*v - 2)
    return w

def get_quartic_pts(u: Rational, max_pt: int=100_00_000, D2=None) -> None|list:
    """Get points on quartic for u or on given quartic.
    """
    if D2 is None:
        D2 = u_to_quartic(u)

    # If both D2 and -D2 return nothing too quickly, return None
    run_l = run_m = False
    d = 100_000
    n = 5 * d
    time0e = datetime.now()
    p = pari(D2).hyperellratpoints([n, [4, d]], 0) 
    time1e = datetime.now()
    l1 = list(p)
    #print(f'{u}, l1 {len(l1)}, {time1e-time0e}')
    if time1e-time0e >= timedelta(milliseconds=1) or 0 < len(l1): 
        run_l = True

    time0e = datetime.now()
    p = pari(-D2).hyperellratpoints([n, [4, d]], 0) 
    time1e = datetime.now()
    m1 = list(p)
    #print(f'{u}, l1m {len(l1m)}, {time1e-time0e}')
    if time1e-time0e >= timedelta(milliseconds=1) or 0 < len(m1): 
        run_m = True
    if not run_l and not run_m: return None

    # Search for more pts
    print(f'Searching for pts on {u} quartic')

    l2 = l3 = l4 = l5 = l6 = m2 = m3 = m4 = m5 = m6 = []
    if max_pt > d:
        s = d + 2
        d = min(1_000_000, max_pt)
        n = 3 * d
        if run_l:
            p = pari(D2).hyperellratpoints([n, [s, d]], 0)
            l2 = list(p)
        if run_m:
            p = pari(-D2).hyperellratpoints([n, [s, d]], 0)
            m2 = list(p)
        if max_pt > d:
            s = d + 2
            d = min(10_000_000, max_pt)
            n = 2 * d
            if run_l:
                p = pari(D2).hyperellratpoints([n, [s, d]], 0)
                l3 = list(p)
            if run_m:
                p = pari(-D2).hyperellratpoints([n, [s, d]], 0)
                m3 = list(p)
            if max_pt > d:
                s = d + 2
                d = min(50_000_000, max_pt)
                n = d
                if run_l:
                    p = pari(D2).hyperellratpoints([n, [s, d]], 0)
                    l4 = list(p)
                if run_m:
                    p = pari(-D2).hyperellratpoints([n, [s, d]], 0)
                    m4 = list(p)
                if max_pt > d:
                    s = d + 2
                    d = min(80_000_000, max_pt)
                    n = d
                    if run_l:
                        p = pari(D2).hyperellratpoints([n, [s, d]], 0)
                        l5 = list(p)
                    if run_m:
                        p = pari(-D2).hyperellratpoints([n, [s, d]], 0)
                        m5 = list(p)
                    if max_pt > d:
                        s = d + 2
                        d = max(100_000_000, max_pt)
                        n = d
                        if run_l:
                            p = pari(D2).hyperellratpoints([n, [s, d]], 0)
                            l6 = list(p)
                        if run_m:
                            p = pari(-D2).hyperellratpoints([n, [s, d]], 0)
                            m6 = list(p)

    l_pts = l1 + l2 + l3 + l4 + l5 + l6
    m_pts = m1 + m2 + m3 + m4 + m5 + m6
    pts = [(QQ(v), QQ(D)) for (v, D) in l_pts + m_pts]
    return pts

def find_uvD_pts(first_ud: int, last_ud: int, first_un: int, last_un: int,
        max_pt: int=int(1e8), max_d: int=int(1e27)) -> None | tuple:
    """ 
    Search for solution pairs using hyperellratpoints with:
        ud in range(4<=first_ud, last_ud+1, 4);
        un in range(1<=first_un, last_un+1, 2);
        (u_num, u_den) in ((un, ud), (-un, ud), (ud, un), (ud, -un));
        v a rational point on the quartic defined by u up to max_pt
    Check u_den / u_num for quadratic obstruction
        to the solutions of the quadratics for y^2 and t^2.
    """
    # Check input parameters for sanity.
    assert first_ud % 4 == 0
    assert first_un % 2 == 1
    assert 4 <= first_ud <= last_ud
    assert 1 <= first_un <= last_un
    assert 1 <= max_pt

    # Look up index of known solution with denominator.
    d_to_known_inx = {val['abcd'][-1]: inx 
            for inx, val in enumerate(known.values(), start=1)}

    # Declare search counts
    u_hits = D_hits = big_hits = known_hits = 0
    found_knowns = set() # indexes of known solutions found in search
    found_new = False

    # Search for points and count results.
    start = datetime.now()
    for ud in range(first_ud, last_ud + 1, 4):
        for un in range(first_un, last_un + 1, 2):
            if gcd(un, ud) != 1: continue
            uQ = QQ(un) / ud
            for u in (uQ, -uQ, uQ.inverse(), -uQ.inverse()):
                # Quadratics for y^2 and t^2 on inverse must be solvable
                if not check_yt(u.inverse()): continue
                u_hits += 1
                pts = get_quartic_pts(u, max_pt)
                if pts is None: continue
                if 0 == len(pts): continue
                for v, D in pts:
                    if D < 0: continue
                    D_hits += 1
                    w = uv_to_w(u, v)
                    print(f'u {u}, v {v} -> D {D}, w {w}', flush=True)
                    pair = uvD_to_xyz(u, v, D)
                    for xyz in pair:
                        d = lcm([x.denom() for x in xyz])
                        a, b, c = abc = sorted(abs(x) * d for x in xyz)
                        assert a**4 + b**4 + c**4 == d**4
                        if d >= max_d:
                            big_hits += 1
                            print(f'\tbig {float(d):.4e}')
                            continue
                        if d in d_to_known_inx:
                            known_hits += 1
                            inx = d_to_known_inx[d]
                            found_knowns.add(inx)
                            print(f'\tknown #{inx}: {d}', flush=True)
                            continue
                        print(f'\n\nNEW {u}, {v} -> {d}; {c}, {b}, {a}',
                              flush=True)
                        found_new = True
                        break
                    if found_new: break
                if found_new: break
            if found_new: break
        if found_new: break

    # Report results and return if found new solution.
    elapsed = datetime.now() - start
    print(f'hits: u {u_hits}, D {D_hits}, big {big_hits}, known {known_hits}'
        f'\nknowns {len(found_knowns)}: {sorted(found_knowns)}')
    print(f'elapsed: {elapsed}')
    if found_new:
        return u, v, xyz
"""
find_uvD_pts(4, 4, 201, 201, 1010)
1: u 201/4
201/4, -136/133 -> D 1416600375/141512
	known #12: (27450160/156646737, 146627384/156646737, 108644015/156646737)
	known #5: (-12552200/16003017, -4479031/16003017, -14173720/16003017)
201/4, -1005/568 -> D 87999215295/5161984
	known #5: (4479031/16003017, 12552200/16003017, 14173720/16003017)
	known #12: (-146627384/156646737, -27450160/156646737, -108644015/156646737)
2: u -201/4
	no points for u -201/4
hits: u 2, D 2, big 0, known 4
knowns 2: [5, 12]
elapsed: 0.1065s

find_uvD_pts(200, 200, 1, 2001, 100_000)
-1617/200, -5/8 -> D 708019737/2560000
	known #13: (186668000/589845921, 582665296/589845921, 260052385/589845921)
	known #14: (-219076465/638523249, -275156240/638523249, -630662624/638523249)
-1617/200, -34272/4885 -> D 2365840162968393/477264500000
	known #14: (275156240/638523249, 219076465/638523249, 630662624/638523249)
	known #13: (-582665296/589845921, -186668000/589845921, -260052385/589845921)
1873/200, -51416/9425 -> D 4964755565640087/1776612500000
	known #41: (487814048600/26969608212297, 8528631804200/26969608212297, 26901926181047/26969608212297)
	known #52: (-103028409596553328/103117303193818953, 24975412054750025/103117303193818953, -4092004076331400/103117303193818953)
hits: u 149, D 3, big 0, known 6
knowns 4: [13, 14, 41, 52]
elapsed: 22.7084
"""


def DEBUG(*args):
    import pdb; pdb.set_trace()
    pass; pass; pass # opportumity to debug
"""
python -m elliptical.solutions_curves DEBUG 'QQ(4/6)' 4/6
(Pdb) args
args = (2/3, 2/3)
(Pdb) 4/6
0.6666666666666666
"""

if __name__ == "__main__":
    import sys
    command = eval(sys.argv[1])
    print(command)
    args = list(map(sage_eval, sys.argv[2:]))
    result = command(*args)
    print(result)

"""
TODO
How do EC from pairs compare to EC from solutions_curves.log?
    python -m elliptical.solutions_curves get_optimized_rational_points -20/9
    q_res = make_quartic(QQ(20/-9) with different quad_xy
    get_quartic_pts(QQ(20/-9), 1_000_000, q_res[2])
        with a dozen different quad_xy
    e_res = model_quartic_as_elliptic_curve(q_res[2], k0)
    factor(C0.conductor())
vs  JD2 = u_to_D_to_EC
The quartics differ, but the EC remain the same.

Find places in solutions_curves with couldn't find quartic points.
And use get_quartic_pts(QQ(20/-9), 1_000_000, q_res[2])

How go from pnt on EC back to pnt on quartic and solution to abcd?

Methods of building EC
Use pieza_in_pairs.py
    u_to_D_to_EC(u)

Use solutions_curves.py
    make_quartic(mn: Rational, quad_xy: tuple)

Use curves_tomita.py
    model_quartic_as_elliptic_curve(quartic_poly, quartic_x: Rational)
Use modular.log
# The Sage-native way to convert a quartic to an elliptic curve
R.<x,y> = QQ[]
f = x^4 + x + 1  # your quartic
C = Curve(y^2 - f)
E = Jacobian(C)
print(E.minimal_model())

-------------

>>> abcd_to_small_u((5_870_000, 8_282_543, 11_289_040, 12_197_457), 9e99)
[-400/37, -1010819791893/158785763120, -338000022077/104832225800, -2433/920, 
-35798568240/28609248113, -93/80, -84237/359800, 20632147/117135680, 
8685847/22963880, 450668400/124346123, 11502160/2925527, 1867333/457280]
>>> JD2 = u_to_D_to_EC(QQ(-93/80))
>>> JD2[0]
Elliptic Curve defined by y^2 = x^3 + 1687750917943790881*x - 294299265667029867546450078 over Rational Field
>>> JD2[1]
876967441/40960000*v^4 - 1160148721/10240000*v^3 + 5260917/32000*v^2 - 930349361/5120000*v - 540248559/10240000
>>> JD2[0].conductor()
31173291851033505900611518136
>>> from sage.all import factor
>>> factor(JD2[0].conductor())
2^3 * 17 * 41 * 73 * 89 * 241 * 1249 * 2137 * 2887 * 6569 * 70537
>>> JD2[0].torsion_points()
[(0 : 1 : 0), (171390639 : 0 : 1)]
>>> JD2[0].hyperelliptic_polynomials()[0].roots()[0][0]
171390639
>>> pari(JD2[0]).ellrank(15)
[2, 2, 2, [[18821840650969/8649, 94094386506723420026/804357], [2216934161202292017321/77272324441, 104489172057583918416101354132502/21480083475784739]]]
ellrank returns different values each time.
[2, 2, 2, [[1057203921, 51687652673778], [7086871502336635716361489/30329632086478209, 56104185538669647074758265030211769294/5282028171881234074661823]]]
 
>> g = JD2[0].gens(algorithm="pari", pari_effort=15)
[(7086871502336635716361489/30329632086478209 : 56104185538669647074758265030211769294/5282028171881234074661823 : 1), (18821840650969/8649 : 94094386506723420026/804357 : 1)]
but sage cache the values.

[r,R,s,V] = ellrank(E)
W = ellsaturation(E, V, 100)
takes too long when V is longer than the rank.

How to do LLL basis reduction on the generators to get smaller height points.
# 1. Get the height matrix from your rank-2 (or rank-k) curve
M = E.height_pairing_matrix(gens)

# 2. Call PARI's qflll on that matrix
# We convert the Sage matrix to a PARI object first
pari_M = pari(M)
T = pari_M.qflllgram() 

# 3. Apply the transformation T to your generators
# If gens = [P1, P2], and T = [[t11, t12], [t21, t22]]
# New P1 = t11*P1 + t21*P2, etc.
new_gens = []
for col in range(len(gens)):
    # Combine original generators using the coefficients in the T matrix
    new_P = sum(int(T[i, col]) * gens[i] for i in range(len(gens)))
    new_gens.append(new_P)

# 4. new_gens now contains the LLL-reduced basis (shortest heights)

"""