"""
Code to explore how to manipulate elliptic curves like Seiji Tomita does.
See my notes in curves_tomita.md
and Tomita's notes in http://www.maroon.dti.ne.jp/fermat/dioph4e.html
"""
from sage.all import EllipticCurve
from elliptical.solutions_curves import *

def abcd_to_xyz(abcd: tuple) -> tuple:
    """Convert (a, b, c, d) representing a^4 + b^4 + c^4 = d^4 
    to rational form (x, y, z) representing x^4 + y^4 + z^4 = 1
    """
    A, B, C, D = map(QQ, abcd)
    x, y, z = A/D, B/D, C/D
    return x, y, z

def xyz_to_u(x :Rational, y :Rational, z :Rational) -> Rational:
    """Convert (x, y, z) representing x^4 + y^4 + z^4 = 1
    to single parameter u.
    """
    x_y = x - y
    u = (x_y*x_y - z*z - 1) / (x*x - x*y + y*y + x_y)
    return u

def abcd_to_u_trace(abcd: tuple) -> set:
    """Convert (a, b, c, d) to (x, y, z) and then consider
    3 permutations and 8 combinations of signs on (x, y, z).
    Since the u, v, w formulas rotate (x, y, z), only swap y, z.
    Trace which u come from permutations. Return sorted set.
    """
    u_set = set()
    xx, yy, zz = abcd_to_xyz(abcd)
    for x, y, z in ((xx, yy, zz), (zz, xx, yy), (zz, yy, xx),):
        for sx in (x, -x):
            for sy in (y, -y):
                for sz in (z,):
                    u = xyz_to_u(sx, sy, sz)
                    u_set.add(u)
                    print(f'u {u} from {sx, sy, sz}')
    return set(sorted(u_set))


def model_quartic_as_elliptic_curve(quartic_poly, quartic_x: Rational):
    """Convert parameterized quartic to Elliptic Curve
    """
    y_sq = quartic_poly.subs(k=quartic_x)
    if not y_sq.is_square():
        quartic_poly = -quartic_poly
        y_sq = quartic_poly.subs(k=quartic_x)        
    y0 = y_sq.sqrt()
   
    # u = k - k0 is point shifted to origin
    R = PolynomialRing(QQ, 'u')
    u = R.gen()
    shifted_poly = quartic_poly.subs({quartic_poly.parent().gen(): u + quartic_x})
    coeffs = shifted_poly.list()
    
    # Standard form: Y^2 = a*u^4 + b*u^3 + c*u^2 + d*u + e^2
    # Ensure we have 5 coefficients
    while len(coeffs) < 5: coeffs.append(0)
    e2, d, c, b, a = coeffs
    #a, b, c, d, e2 = quartic_poly.list()

    #print('Elliptic coeffs with e2\n', coeffs)
    e = y0 # The sqrt of the constant term
    
    # 2. Coefficients of the Weierstrass form y^2 + a1*xy + a3*y = x^3 + a2*x^2 + a4*x + a6
    # These formulas come from the Mordell/Cassels transformation
    a1 = d/e
    a2 = c - (d**2)/(4*e**2)
    a3 = 2*e*b
    a4 = -4*e**2 * a
    a6 = a2 * a4
    
    E_orig = EllipticCurve([a1, a2, a3, a4, a6])
    E_jmin = EllipticCurve(j=E_orig.j_invariant()) # often shorter, but quadratic twist
    if E_jmin.is_isomorphic(E_orig): E_min = E_jmin
    else: 
        E_min = E_orig.short_weierstrass_model().minimal_model()
        assert E_min.is_isomorphic(E_orig)
    return E_min, E_orig

def elliptic_point_to_k(P, E_res: tuple, 
            quartic_poly: Polynomial_rational_flint, k0: Rational) -> None|Rational:
    """
    Map a (new) Elliptic Curve Point back for parameter k on quartic_poly
    Handle sign flips and quadratics twists.
    """
    # Calcualate y for quartic k0
    poly = quartic_poly
    y_sq = poly.subs({poly.parent().gen(): k0})
    if not y_sq.is_square():
        poly = -poly
        y_sq = poly.subs({poly.parent().gen(): k0})
        if not y_sq.is_square():
            import pdb; pdb.set_trace()
            pass # neither sign poly is square
    y = y_sq.sqrt()

    # Map EC point back to Weierstrass to u
    E_min, E_orig = E_res
    assert E_orig.is_isomorphic(E_min)
    iso = E_min.isomorphism_to(E_orig)
    X, Y, Z = P_w = iso(P)
    if P_w.is_zero(): return k0
    
    if y == 0:
        k = (1 / (X/Z)) + k0 # k = 1/x
    elif Y == 0: 
        k = k0 # map root of cubic, back to quartic root
    else:
        R_u = PolynomialRing(QQ, 'u')
        u_var = R_u.gen()
        shifted_poly = poly.subs({quartic_poly.parent().gen(): u_var + k0})
        coeffs = shifted_poly.list()
        while len(coeffs) < 5: coeffs.append(0)
        y2, d, c, b, a = coeffs
        a2 = c - (d**2)/(4*y_sq) 
        u_val = (2 * y * (X/Z + a2)) / (Y/Z)
        k = u_val + k0
    
    if not poly(k).is_square():
        import pdb; pdb.set_trace()
        pass # new quartic k is not square
    return k

def search_small_EC_points(e_res: tuple, quartic_poly: Polynomial_rational_flint, 
                           k0: Rational, mn: Rational, 
                           quad_xy: tuple):
    # Check that EC root point maps back to quartic
    E_min = e_res[0]
    p0 = E_min(E_min.hyperelliptic_polynomials()[0].roots()[0][0], 0)
    k = elliptic_point_to_k(p0, e_res, quartic_poly, k0)
    print(f'p0 = {p0}, k = {k}')
    abc, d = k_to_abcd(mn, quad_xy, k) # Check that p0 works
    print(f'\td {d}')

    # Check that generators map back to quartic
    E = pari(e_res[0])
    g2 = set() # gens in (x, y) form
    for effort in (2, 2, 2, 2, 3, 3, 3, 4, 4, 5):
        rank = list(E.ellrank(effort))
        g2.update(rank[-1])
        if rank[0] != rank[1]: continue
        print(f' effort {effort}, rank {rank}')
        break
    print(f'Found {len(g2)} gens')
    gs = {E_min((g[0], g[1], 1)) for g in g2} # convert to (x, y, z) form
    for g in gs:
        k = elliptic_point_to_k(g, e_res, quartic_poly, k0)
        print(f'g = {g}, k = {k}')
        abc, d = k_to_abcd(mn, quad_xy, k) # Check that g works
        print(f'\td {d}')
    
    # Augment generators
    gs |= {-g for g in gs}
    gs |= {2*g for g in gs}
    gs |= {g1 + g2 for g1 in gs for g2 in gs}
    gs |= {g1 + g2 + g3 for g1 in gs for g2 in gs for g3 in gs}

    # Use p0 and augmented generators to build list of points
    p_set = set((p0,)) | gs
    p_set |= {p0 + g for g in gs}
    p_list = sorted(p_set)
    print(f'Using {len(p_list)} points on EC from {len(gs)} generators')

    # Generate solutions from points
    d_count = 0
    d_ok_set = set()
    d_big_set = set()
    for p in p_list:
        k = elliptic_point_to_k(p, e_res, quartic_poly, k0)
        abc, d = k_to_abcd(mn, quad_xy, k)
        d_count += 1
        if d < 1e27:
            if d in d_ok_set: continue
            d_ok_set.add(d)
            print(f'd={d}, abc={abc}, EC p {p} -> k {k}')
        else:
            d_big_set.add(d)
    print(f'{d_count} d, {len(d_ok_set)} ok, {len(d_big_set)} big')
    return sorted(d_ok_set)

def DEBUG():
    import pdb; pdb.set_trace()
    pass; pass; pass # opportumity to debug

""" TODO
"""

if __name__ == "__main__":
    import sys
    command = eval(sys.argv[1])
    print(command)
    args = list(map(int, sys.argv[2:]))
    result = command(*args)
    print(result)

