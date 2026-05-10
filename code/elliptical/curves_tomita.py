"""
Code to explore how to manipulate elliptic curves like Seiji Tomita does.
See my notes in curves_tomita.md
and Tomita's notes in http://www.maroon.dti.ne.jp/fermat/dioph4e.html
"""
import itertools
from elliptical.solutions_curves import *
from sage.all import EllipticCurve

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
                           k0: Rational, mn: Rational, quad_xy: tuple,
                           coeff_lim: int=3, k_h_lim: int=int(1e10), 
                           verbose=False, max_d: int=int(1e27)) -> list:
    """ Extract elliptic curve. Walk the EC for small points, 
    and map back to quart points k limied by height.
    """
    # Check that EC root point maps back to quartic
    E = e_res[0]
    k_set = set()

    # Check that torsion points map back to quartic
    tor = E.torsion_points()
    for t in tor:
        k = elliptic_point_to_k(t, e_res, quartic_poly, k0)
        k_set.add(k)
        if verbose: print(f't = {t}, k = {k}')
        abc, d = k_to_abcd(mn, quad_xy, k) # Check that t works
        if verbose: print(f'\td {d}')

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
        k = elliptic_point_to_k(g, e_res, quartic_poly, k0)
        k_set.add(k)
        if verbose: print(f'g = {g}, k = {k}')
        abc, d = k_to_abcd(mn, quad_xy, k) # Check that g works
        if verbose: print(f'\td {d}')

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

    # Map back to k from points
    for p in sorted(p_set):
        k = elliptic_point_to_k(p, e_res, quartic_poly, k0)
        k_set.add(k)
    
    # Restrict k by height
    ks = [(k.height(), k) for k in k_set if k.height() < k_h_lim]
    if verbose: print(f'Restrict to {len(ks)} k with height < {k_h_lim}')
    ks = [k for h, k in sorted(ks)]

    # Report abcd solutions for points. Index known with denominator
    d_to_known_inx = {val['abcd'][-1]: inx 
        for inx, val in enumerate(known.values(), start=1)}

    # Generate solutions from points
    big_hits = known_hits = 0
    found_knowns = set()
    for k in ks:
        abc, d = k_to_abcd(mn, quad_xy, k)
        if d >= max_d:
            big_hits += 1
            if verbose: print(f'\tk {k} -> big d {float(d):.5e}')
            continue
        if d in d_to_known_inx:
            known_hits += 1
            inx = d_to_known_inx[d]
            found_knowns.add(inx)
            if verbose: print(f'\tk {k} -> known #{inx}: {d}')
            continue
        a, b, c = abc
        print(f'\n\nNEW: k {k} -> {d}; {c}, {b}, {a}', flush=True)
        return k, d, abc
    found_knowns = sorted(found_knowns)
    print(f'big {big_hits}, known {known_hits}'
        f'\nknowns {len(found_knowns)}: {found_knowns}')
    return found_knowns

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

