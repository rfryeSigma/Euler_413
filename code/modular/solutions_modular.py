"""
Explore whether modular curves can help find new solutions to a^4 + b^4 + c^4 = d^4.
"""
import csv
from datetime import datetime, timedelta
import itertools
from math import gcd
from pdb import set_trace, runcall
from solutions import known
from sage.all import help, oo, pari, sage_eval, \
    QQ, RR, Rational, hilbert_symbol, lcm, \
    DiagonalQuadraticForm, EllipticCurve, PolynomialRing
from sage.rings.polynomial.polynomial_rational_flint import Polynomial_rational_flint
from sage.modules.free_module_element import vector
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

def u_to_E_min_map(u, v0):
    """ Generate D2 quartic, minimal elliptic curve and
    map from elliptical points to quartic points.
    """
    D2 = u_to_quartic(u) # make quartic polynomial in v

    # Shift the quartic so v0 is at the origin: u = v - v0
    R = PolynomialRing(QQ, 'u')
    u = R.gen()
    shifted_poly = D2.subs({D2.parent().gen(): u + v0})
    coeffs = shifted_poly.list()
    while len(coeffs) < 5: coeffs.append(0)
    
    # y^2 = a*u^4 + b*u^3 + c*u^2 + d*u + e^2
    e2, d, c, b, a = coeffs
    e = e2.sqrt()
    
    # Map to Weierstrass: y^2 + a1*xy + a3*y = x^3 + a2*x^2 + a4*x + a6
    a1 = d/e
    a2 = c - (d**2)/(4*e**2)
    a3 = 2*e*b
    a4 = -4*e**2 * a
    a6 = a2 * a4
    
    E_orig = EllipticCurve([a1, a2, a3, a4, a6])
    E_min = E_orig.minimal_model()
    iso = E_min.isomorphism_to(E_orig)
    iso_inv = E_orig.isomorphism_to(E_min)
    
    # Return everything needed to bridge back
    bridge = {
        'v0': v0,
        'D2': D2,
        'e': e,
        'a2': a2,
        'iso_to_orig': iso,
        'iso_from_orig': iso_inv,
        'shifted_poly': shifted_poly
    }
    return E_min, bridge

def map_point_to_v(P, bridge):
    """ Map a (new) Elliptic Curve Point back for parameter k on D2
    Inverts the transformation: E_min -> E_orig -> u -> v
    """
    # Map from Minimal Model back to the specific Weierstrass model
    e = bridge['e']
    a2 = bridge['a2']
    v0 = bridge['v0']
    if P.is_zero(): return v0
    P_w = bridge['iso_to_orig'](P)
    if P_w.is_zero(): return v0

    X, Y = P_w[0], P_w[1]
    if Y == 0: return v0 # Maps to the seed point (or a root)

    # Apply the inverse Mordell transformation     
    u_val = (2 * e * (X + a2)) / Y
    return u_val + v0

def map_v_to_point(v, bridge):
    """
    Map a rational value v from the quartic D2 back to a point on the EC.
    Returns a point on E_min.
    """
    # Shift v to the origin: u = v - v0
    v0 = bridge['v0']
    u = v - v0
    
    # The seed point v0 always maps to the Point at Infinity
    iso_inv = bridge['iso_from_orig']
    if u == 0:
        # iso_inv maps the identity of E_orig to the identity of E_min
        return iso_inv.domain()(0) 
    
    # Determine y such that y^2 = D2(v)
    y_sq = bridge['D2'](v)
    if not y_sq.is_square():
        # If D2 was negated during bridge creation, account for that
        y_sq = -y_sq
        assert y_sq.is_square()
    y = y_sq.sqrt()
    
    # Extract coefficients from the shifted poly for the map
    # shifted_poly: y^2 = a*u^4 + b*u^3 + c*u^2 + d*u + e^2
    coeffs = bridge['shifted_poly'].list()
    while len(coeffs) < 5: coeffs.append(0)
    e2, d, c, b, a = coeffs
    
    # Apply the forward Mordell transformation to E_orig
    # These formulas map (u, y) on the quartic to (X, Y) on the Weierstrass form
    e = bridge['e']
    X = (2*e*(y + e) + d*u) / (u**2)
    Y = (4*e**2*(y + e) + 2*e*(d*u + c*u**2) - (d**2 * u**2 / (2*e))) / (u**3)
    
    # Create the point on E_orig
    E_orig = iso_inv.domain()
    P_orig = E_orig(X, Y)
    
    # 6. Map to E_min
    return iso_inv(P_orig)

def decompose_point(E, P, gens, torsion_points):
    """
    Finds n_i such that P = sum(n_i * G_i) + T.
    Uses .height() which is compatible with both Rational and Number Field points.
    """
    from sage.matrix.constructor import matrix

    # Force P onto the exact curve object E
    P = E(P) 

    for T in torsion_points:
        P_test = P - T
        if P_test.is_zero():
            return [0] * len(gens), T
            
        # Build the Gram Matrix using the Bilinear Pairing formula
        # <P, Q> = (h(P+Q) - h(P) - h(Q)) / 2
        dim = len(gens)
        M = matrix(RR, dim, dim)
        for i in range(dim):
            for j in range(dim):
                h_gi = (gens[i]).height()
                h_gj = (gens[j]).height()
                h_sum = (gens[i] + gens[j]).height()
                M[i,j] = (h_sum - h_gi - h_gj) / 2

        # Build the B vector for P_test
        h_ptest = P_test.height()
        B_list = []
        for G in gens:
            h_g = G.height()
            h_comb = (P_test + G).height()
            B_list.append((h_comb - h_ptest - h_g) / 2)
        B = vector(RR, B_list)
        
        try:
            # Solve M * n = B
            coeffs = M.solve_right(B)
            n_vals = [round(float(c)) for c in coeffs]
            
            # 4. Final Verification
            verification = sum(n_vals[i] * gens[i] for i in range(len(gens)))
            if verification == P_test:
                return n_vals, T
        except Exception:
            continue   
    return None, None

def find_uv_by_EC(u: Rational, v0: Rational, coeff_lim: int=3,
            v_h_lim=1e10, verbose=False):
    """ Build D2 elliptic curve, and return map.
    Walk the EC for small points, and map back to quart points v
    limied by height.
    """
    E, bridge = u_to_E_min_map(u, v0)
    v_set = set()

    # Check that torsion points map back to quartic
    tor = E.torsion_points()
    for t in tor:
        v = map_point_to_v(t, bridge)
        v_set.add(v)
        if verbose: print(f't = {t}, v = {v}')

    # Collect generators
    g2 = set() # gens in (x, y) form
    R = pari(E).ellrankinit() # speed up ellrank
    for effort in (3, 3, 4, 5):
        rank = R.ellrank(2) # initialize saved points
        rank = R.ellrank(effort-1, rank[3]) # don't redo saved points
        rank = R.ellrank(effort, rank[3])
        for g in rank[3]: g2.add(g)
        if verbose: print(f' effort {effort}, rank {rank[:3]}, {len(g2)} gens')
    gs = [E((g[0], g[1], 1)) for g in g2] # convert to (x, y, z) form

    # Check whether EC root point is indepenent of generators
    p0 = E(E.hyperelliptic_polynomials()[0].roots()[0][0], 0)
    if p0 in tor:
        if verbose: print(f'Root {p0} is a torsion point')
    elif E.is_independent([p0] + gs):
        if verbose: print(f'Adding root {p0} to gens')
        gs.append(p0)

    # Guarantee full basis of generators
    gs = E.saturation(gs)[0]
    if verbose: print(f'Found {len(gs)} basis gens')
 
    # Check that generators map back to quartic
    for g in gs:
        v = map_point_to_v(g, bridge)
        v_set.add(v)
        if verbose: print(f'g = {g}, v = {v}')

    # Storage for points: n1*G1 + n2*G2 + n3*G3 + n4*G4
    # With 4 generators and limit 5, this is 11^4 = 14,641 combinations
    ps = [None] * len(tor) * (2 * coeff_lim + 1) ** len(gs)
    p_inx = 0
    coeff_ranges = [range(-coeff_lim, coeff_lim + 1) for _ in range(len(gs))]
    for coeffs in itertools.product(*coeff_ranges):
        # Compute the linear combination n1*G1 + n2*G2 + ...
        P = E(0)
        for i, n in enumerate(coeffs):
            P += n * gs[i]
        
        # Add the torsion points to this combination
        for t in tor:
            ps[p_inx] = t + P
            p_inx += 1
    assert p_inx == len(ps)
    p_set = set(ps)
    if verbose: print(f'Found {len(p_set)} unique points from {len(ps)}')

    # Map back to v from points
    for p in sorted(p_set):
        v = map_point_to_v(p, bridge)
        v_set.update(v)
    
    # Restrict v by height
    vs = [(v.height(), v) for v in v_set if v.height() < v_h_lim]
    return [v for h, v in sorted(vs)]

def test_coeffs(n_gs: int=3, coeff_lim: int=3):
    """ Show cominations of generator indexes and multipliers
    There are (2 * coeff_lim + 1)**n_gs combinations
    but one combo with n = 0 on all generators is unused.
    """
    # With 1 generators and limit 1, this is 3 - 1 combos
    # With 2 generators and limit 1, this is 9 - 1 combos
    # With 1 generators and limit 2, this is 5 - 1 combos
    #
    # With 4 generators and limit 5, this is 11^4 = 14,641 -1 combinations
    coeff_ranges = [range(-coeff_lim, coeff_lim + 1) for _ in range(n_gs)]
    for coeffs in itertools.product(*coeff_ranges):
        print(list(coeffs))
        for i, n in enumerate(coeffs):
            print('\t', i, n)
""" Yes, it generates all combinations of generators and multipliers
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

def solve_v_list(u: Rational, v_list: list, max_d: int=int(1e27)):
    """ Report abcd solutions for u, v pairs
    """
    # Look up index of known solution with denominator.
    d_to_known_inx = {val['abcd'][-1]: inx 
            for inx, val in enumerate(known.values(), start=1)}
    D_hits = big_hits = known_hits = 0
    found_knowns = set() # indexes of known solutions found in search
    found_new = False
    D2 = u_to_quartic(u)
    for v in v_list:
        d2 = D2(v)
        assert d2.is_square()
        D = d2.sqrt()
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
    print(f'big {big_hits}, known {known_hits}'
        f'\nknowns {len(found_knowns)}: {sorted(found_knowns)}')
    if found_new:
        return u, v, xyz

def uv_to_w(u: Rational, v: Rational) -> Rational:
    """The u, v, w derived from x^4 + y^4 + z^4 = 1
    have relationship 2(u+v+w)-uvw-4=0
    """
    w = 2 * (u + v - 2) / (u*v - 2)
    return w

def get_quartic_pts(u: Rational, max_pt: int=100_00_000, D2=None,
    d_list: tuple=(100_000, 1_000_000, 10_000_000, 50_000_000, 80_000_000, 100_000_000),
        n_mult: tuple= (5, 3, 2, 1, 1, 1)) -> None|list:
    """Get points on quartic for u or on given quartic.
    """
    if D2 is None:
        D2 = u_to_quartic(u)

    # If both D2 and -D2 return nothing too quickly, return None
    run_l = run_m = False
    d = d_list[0]
    n = n_mult[0] * d
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
        d = min(d_list[1], max_pt)
        n = n_mult[1] * d
        if run_l:
            p = pari(D2).hyperellratpoints([n, [s, d]], 0)
            l2 = list(p)
        if run_m:
            p = pari(-D2).hyperellratpoints([n, [s, d]], 0)
            m2 = list(p)
        if max_pt > d:
            s = d + 2
            d = min(d_list[2], max_pt)
            n = n_mult[2] * d
            if run_l:
                p = pari(D2).hyperellratpoints([n, [s, d]], 0)
                l3 = list(p)
            if run_m:
                p = pari(-D2).hyperellratpoints([n, [s, d]], 0)
                m3 = list(p)
            if max_pt > d:
                s = d + 2
                d = min(d_list[3], max_pt)
                n = n_mult[3] * d
                if run_l:
                    p = pari(D2).hyperellratpoints([n, [s, d]], 0)
                    l4 = list(p)
                if run_m:
                    p = pari(-D2).hyperellratpoints([n, [s, d]], 0)
                    m4 = list(p)
                if max_pt > d:
                    s = d + 2
                    d = min(d_list[4], max_pt)
                    n = n_mult[4] * d
                    if run_l:
                        p = pari(D2).hyperellratpoints([n, [s, d]], 0)
                        l5 = list(p)
                    if run_m:
                        p = pari(-D2).hyperellratpoints([n, [s, d]], 0)
                        m5 = list(p)
                    if max_pt > d:
                        s = d + 2
                        d = max(d_list[5], max_pt)
                        n = n_mult[5] * d
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
"""

# unfinished
def scan_known(n_known: int=8, max_h: int=100_000_000, verbose=False,
               max_d: int=int(1e27)) -> None:
    """For all abcd in first n_known solutions, generate 12 u.
    For each u with height less than max_h:
        map u -> d, E a4, a6, E conductor, paired v
            Walk E with v for solutions
    ...
    """
    u_info = dict() # dict of info on bounded u
        # ds: list of d from abcd
        # a_46: Elliptic Curve (a4, a6) coefficients
        # conductor: EC conductor
        # vs: v found in EC
    for val in list(known.values())[:n_known]:
        a, b, c, d = abcd = val['abcd']
        if verbose: print(d)
        hu_s = abcd_to_h_u(abcd)
        for h, u in hu_s:
            if h > max_h: continue
            #print(u, check_elkies_rules(u))
            assert check_yt(u.inverse()), f'No y and t for u {u}'
            if verbose: print('\t', u)
            info = u_info[u] = u_info.get(u, dict())
            ds = info['ds'] = info.get('d', list())
            ds.append(d)
            vs = info['vs'] = info.get('v', list())
            D2 = u_to_quartic(u)
            for _, v in hu_s:
                if v == u: continue
                d2 = D2(v)
                if not d2.is_square(): continue
                vs.append(v)
            v0 = vs[0]
            E_min, bridge = u_to_E_min_map(u, v0)
            info['a_46'] = (E_min.a4(), E_min.a6())
            info['conductor'] = E_min.conductor()
            v_list = find_uv_by_EC(u, v0)
            assert solve_v_list(u, v_list) is None, f'{u}: {info}' 

    #print(f'#d {len({info['ds']) for info in u_info.values()})}')
    print(f'#E {len({info['a_46'] for info in u_info.values()})}')
    print(f'#c {len({info['conductor'] for info in u_info.values()})}')
    pass
    return u_info
    """
    python -um modular.solutions_modular scan_known 200
    <function scan_known at 0x1561b9760>
    #E 158
    #c 158
    158
    finished after 53s

    u_info = scan_known(200)
    for (u,info) in u_info.items(): 
        assert solve_v_list(u, find_uv_by_EC(u, info['vs'][0])) is None, f'{u}: {info}' 

    (Pdb) debug scan_known()
    ((Pdb)) len({info['d'] for info in u_info.values()})
    8

    """



def DEBUG(*args):
    set_trace()
    pass; pass; pass # opportumity to debug
"""
python -um modular.solutions_modular DEBUG 'QQ(2^3/6)' 2^3/6
(Pdb) args
args = (4/3, 4/3)
(Pdb) 2^3/6
*** TypeError: unsupported operand type(s) for ^: 'int' and 'float'
(Pdb) 8/6
1.3333333333333333
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

>>> abcd_to_small_u((5_870_000, 8_282_543, 11_289_040, 12_197_457), 9e99)
[-400/37, -1010819791893/158785763120, -338000022077/104832225800, -2433/920, 
-35798568240/28609248113, -93/80, -84237/359800, 20632147/117135680, 
8685847/22963880, 450668400/124346123, 11502160/2925527, 1867333/457280]
>>> JD2 = u_to_D_to_EC(QQ(-93/80))
>>> JD2[0]
Elliptic Curve defined by y^2 = x^3 + 1687750917943790881*x - 294299265667029867546450078 over Rational Field
>>> JD2[0].conductor()
31173291851033505900611518136
>>> from sage.all import factor
>>> factor(JD2[0].conductor())
2^3 * 17 * 41 * 73 * 89 * 241 * 1249 * 2137 * 2887 * 6569 * 70537

"""