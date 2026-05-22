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
    GF, QQ, RR, Rational, hilbert_symbol, kronecker, lcm, \
    DiagonalQuadraticForm, EllipticCurve, PolynomialRing
from sage.rings.polynomial.polynomial_rational_flint import Polynomial_rational_flint
from sage.modules.free_module_element import vector
from timeit import repeat

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

def check_quadratic(D, delta):
    """ Check the Hilbert symbol on equation X^2 - Dy^2 - delta*z^2 = 0
    at infinite prime (Real Solvability)
    and at 2 and odd primes (Local) Solvability)
    """
    return hilbert_symbol(D, delta, -1) == 1 == hilbert_symbol(D, delta, 2)

def check_yt(mn: Rational) -> bool:
    """ Return whether both y^2 and t^2 are solvable and have a point
    """
    y2_coeffs, t2_coeffs = mn_to_xyt_conics(mn)
    for a0, a, b, c in (t2_coeffs, y2_coeffs):
        D = 4 * a * a0
        delta = b**2 - 4 * a * c
        if not check_quadratic(D, delta):
            return False
        Q = DiagonalQuadraticForm(QQ, [1, -D, -delta])
        try: # find a rational point (X, V, Z)
            point = Q.solve()
        except Exception as e: return False
    return True

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

def check_quartic(q_poly, verbose=False,
        p_list=(3,5,7,11,13,17)) -> int:
    """ Check Kronecker symbol `(x|y)` for local solvability of quartic.
    Return prime that fails, or 0 on success.
    """
    for inx, p in enumerate(p_list):
        K = GF(p)
        try: q_mod = q_poly.change_ring(K)
        except ZeroDivisionError: continue

        # Sum of Kronecker symbols for finite points
        s = sum(kronecker(q_mod(x).lift(), p) for x in K)

        # Points at infinity
        inf = kronecker(q_mod.leading_coefficient().lift(), p)
        val = p + 1 + s + inf
        if 0 == val: return p # no prime residue
    return 0

def check_known_uv(file_name: str='solutions_uv.csv') -> None:
    """Check hypothesis that inverses of a u, v in known pairs
    are always y^2 and t^2 solvable.
    Also check that quartics for u, v are solvable
    """
    with open(file_name, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file)
        header = next(reader)
        for row in reader:
            u_n, u_d, v_n, v_d = [int(row[i]) for i in range(4)]
            u = QQ(u_n)/u_d
            assert check_yt(u.inverse()), f'quad u in {row}'
            q4 = u_to_quartic(u)
            assert 0 == check_quartic(q4), f'quart u in {row}'
            assert 0 == check_quartic(-q4), f'-quart u in {row}'
            v = QQ(v_n)/v_d
            assert check_yt(v.inverse()), f'quad v in {row}'
            q4 = u_to_quartic(v)
            assert 0 == check_quartic(q4), f'quart v in {row}'
            assert 0 == check_quartic(-q4), f'-quart v in {row}'
            #print(row, u, v)
""" Success
"""

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
    #E_min = E_orig.minimal_model()
    E_min = E_orig.minimal_model().short_weierstrass_model()
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
            v_h_lim: int=int(1e15), verbose=False) -> list:
    """ Build D2 elliptic curve. Walk the EC for small points, 
    and map back to quart points v limied by height.
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
    if verbose: print(f'Restrict to {len(vs)} v with height < {v_h_lim}')
    return [v for h, v in sorted(vs)]

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
                print(f'\tbig {float(d):.5e}')
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

def uv_to_w(u: Rational, v: Rational) -> Rational:
    """The u, v, w derived from x^4 + y^4 + z^4 = 1
    have relationship 2(u+v+w)-uvw-4=0
    """
    w = 2 * (u + v - 2) / (u*v - 2)
    return w

# Totally broken. Can never find a rational solution.
def hensel_lift(q_poly, target_prec=8, p_list=(19, 23, 29, 31,)):
    """
    Find a rational point on y^2 = q_poly using a local mod p seed
    and rational reconstruction, without invoking any LLL lattice step.
    """
    from sage.all import Integer, IntegerModRing, rational_reconstruction, factor, QQ
    dq_dx_poly = q_poly.derivative()

    for p in p_list:
        K = GF(p)
        q_mod = q_poly.change_ring(K)
        dq_dx_mod = dq_dx_poly.change_ring(K)
        print(K)
        print(q_mod)
        print(dq_dx_mod)

        # Find a safe, non-singular starting point mod p
        for x_val in K:
            y_sq = q_mod(x_val)
            if y_sq == 0 or not y_sq.is_square(): continue
            
            # Ensure the point is non-singular
            if dq_dx_mod(x_val) == 0: continue
            
            x = Integer(x_val)
            y = Integer(y_sq.sqrt())
            print(f'\nTry lifting {x, y} in {K}')
            modulus = p

            # Lift the local point to p-adic precision
            for step in range(1, target_prec):
                modulus *= p
                R = IntegerModRing(modulus)
                q_R = q_poly.change_ring(R)
                dq_R = dq_dx_poly.change_ring(R)
                xR = R(x)
                yR = R(y)

                F = yR * yR - q_R(xR)
                if F != 0:
                    dFdx = -dq_R(xR)
                    try: 
                        dx = -F * ~dFdx
                    except ZeroDivisionError: 
                        break # Singular path, abandon this root
                    x = Integer(R(x + dx))
                    xR = R(x)

                qx = q_R(xR)
                err = qx - yR * yR
                if err != 0:
                    dy = err * ~(R(2 * y))
                    y = Integer(R(y + dy))
                
                print(f'{step}: (x,y) lifted to {(x, y)} in {factor(modulus)}')

                # --- RECONSTRUCTION ATTEMPT ---
                # Attempt reconstruction when we have gathered enough p-adic precision
                if step >= 5: 
                    try:
                        x_rat = rational_reconstruction(x, modulus)
                        x_cand = QQ(x_rat)
                        q_val = q_poly(x_cand)
                        if q_val >= 0 and q_val.is_square():
                            print(f'>>> Successfully reconstructed rational point!')
                            return x_cand, q_val.sqrt()
                    except (ArithmeticError, ZeroDivisionError):
                        continue # Keep lifting, not enough precision yet
                        
    return None

# --- Pure Python Setup & Execution Example ---
def run_hensel_lift():
    from sage.all import Integers,rational_reconstruction
    # Define the polynomial ring without using Sage syntax shortcuts
    R = PolynomialRing(QQ, 'v')
    v = R.gen()
    
    # Pre-test the reconsruction
    m = 19**8
    u = QQ(-9)/20
    q_poly = u_to_quartic(u)
    v = QQ(-1041)/320
    assert q_poly(v).is_square()
    assert rational_reconstruction(Integers(m)(v), m) == v
    v = QQ(-1425)/412
    assert q_poly(v).is_square()
    assert rational_reconstruction(Integers(m)(v), m) == v

    # Run the solver
    result = hensel_lift(q_poly)
    print(f"Resulting Point: {result}")
"""
>>> run_hensel_lift()
Try lifting (1, 1) in Finite Field of size 7
Try reconstruct {x} mod {modulus}
Try lifting (3, 1) in Finite Field of size 7
Try reconstruct {x} mod {modulus}
Try lifting (4, 3) in Finite Field of size 7
Try lifting (0, 5) in Finite Field of size 11
Try reconstruct {x} mod {modulus}
Try lifting (1, 3) in Finite Field of size 11
...
Try lifting (2, 14) in Finite Field of size 29
...
Try lifting (23, 14) in Finite Field of size 29
Try reconstruct {x} mod {modulus}
Try lifting (26, 4) in Finite Field of size 29
Try reconstruct {x} mod {modulus}
Resulting Point: None

You are entirely right to be skeptical, and your frustration is 100% justified. 
Every time an LLM is asked to write an "Elkies lattice solver" or a "$p$-adic 
Hensel lift point-finder" for a quartic curve, it invariably manufactures 
a hybrid mathematical phantom.AI models tend to blend three entirely 
different concepts: Schoof-Elkies-Atkin (SÉA) (which uses Hensel lifting 
on modular equations over finite fields to count points), Coppersmith’s 
Method (which uses LLL to find small integer roots of polynomials modulo $M$), 
and Noam Elkies' actual 2000 paper on finding rational points near curves. 
The result is always a flawed script that tries to linearly lift a non-linear 
curve from a random local root mod $p$.A standard Hensel lift is a purely 
local tool; it converges to a $p$-adic number in $mathbb{Q}_p$. Because 
the set of $p$-adic points on your curve is uncountable, a blind local 
lift will almost always shoot off into a transcendental $p$-adic space, 
completely missing the discrete, global rational points you are looking for.

"""

def get_quartic_pts(u: Rational, max_pt: int=100_000, D2=None, verbose=True,
        d_list: tuple=(100_000, 1_000_000, 10_000_000, 50_000_000, 
                       80_000_000, 100_000_000),
        n_mult: tuple= (10, 5, 3, 1, 1, 1)) -> None|list:
    """Get points on quartic for u or on given quartic.
    Return None if not solvable
    Return any points found immediately.
    """
    if D2 is None:
        D2 = u_to_quartic(u)

    # If D2 or -D2 solvable, run them, else quit
    run_p = (0 == check_quartic(D2))
    run_m = (0 == check_quartic(-D2))
    if not run_p and not run_m: return None
    flags = '+-' if run_p and run_m else '+' if run_p else '-'
    if verbose: print(f'Searching for {flags} pts on {u} quartic')

    def search_range(r_inx: int, d: int) -> int|list:
        if max_pt <= d: return []
        s = d + 2
        d = min(d_list[r_inx], max_pt)
        n = n_mult[r_inx] * d
        if run_p:
            p = pari(D2).hyperellratpoints([n, [s, d]], 1)
            lp = list(p)
            if 0 < len(lp):
                pts = [(QQ(v), QQ(D)) for (v, D) in lp]
                return pts
        if run_m:
            p = pari(-D2).hyperellratpoints([n, [s, d]], 1)
            lp = list(p)
            if 0 < len(lp):
                pts = [(QQ(v), QQ(D)) for (v, D) in lp]
                return pts
        return d
    
    d = 2 # so first s = 4
    for r_inx in range(len(d_list)):
        res = search_range(r_inx, d)
        if isinstance(res, list): return res
        d = res
    if max_pt > d and verbose: print(f'Search limit {d} < {max_pt}')
    return []

def find_uvD_pts(first_ud: int, last_ud: int, first_un: int, last_un: int,
        max_pt: int=int(1e8), step_u: int=4, max_d: int=int(1e27)) -> None | tuple:
    """ 
    Search for solution pairs using hyperellratpoints with:
        ud in range(4<=first_ud, last_ud+1, step_u);
        un in range(1<=first_un, last_un+1, 2);
        (u_num, u_den) in ((un, ud), (-un, ud), (ud, un), (ud, -un));
        v a rational point on the quartic defined by u up to max_pt
    Check u_den / u_num for quadratic obstruction
        to the solutions of the quadratics for y^2 and t^2.
    """
    # Check input parameters for sanity.
    if 4 == step_u: 
        assert first_ud % 4 == 0
        assert 4 <= first_ud
    assert first_ud <= last_ud
    assert first_un % 2 == 1
    assert 1 <= first_un <= last_un
    assert 1 <= max_pt

    # Look up index of known solutions with denominator.
    d_to_known_inx = {val['abcd'][-1]: inx 
            for inx, val in enumerate(known.values(), start=1)}

    # Declare search counts
    u_hits = D_hits = big_hits = known_hits = 0
    found_knowns = set() # indexes of known solutions found in search
    found_new = False

    # Search for points and count results.
    start = datetime.now()
    for ud in range(first_ud, last_ud + 1, step_u):
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
    print(f'hits: u {u_hits}, D {D_hits}, big {big_hits}, known {known_hits}')
    if 0 < known_hits: 
        print(f'knowns {len(found_knowns)}: {sorted(found_knowns)}')
    print(f'elapsed: {elapsed}')
    if found_new: return u, v, xyz
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

def scan_known_uv(known_s: int=1, known_e: int=9999, max_h: int=100_000_000, 
            verbose=False, max_d: int=int(1e27)) -> dict:
    """For range of known abcd solutions, generate 12 u.
    For each u with height less than max_h:
        map u -> d, E a4, a6, E conductor, paired v
            Walk E with v for solutions
    Report by u: #Elliptic Curves, #Conductors
    Return u_into.
    """
    u_info = dict() # dict of info on bounded u
        # ds: list of d from abcd
        # a_46: Elliptic Curve (a4, a6) coefficients
        # conductor: EC conductor
        # vs: v found in EC
    for val in list(known.values())[known_s - 1: known_e]:
        a, b, c, d = abcd = val['abcd']
        if verbose: print(d)
        hu_s = abcd_to_h_u(abcd)
        for h, u in hu_s:
            if h > max_h: continue
            #print(u, check_elkies_rules(u))
            assert check_yt(u.inverse()), f'No y^2 and t^2 for u {u}'
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
