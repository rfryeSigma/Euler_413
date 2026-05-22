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

def abcd_to_h_mn(abcd: tuple) -> list:
    """Convert (a, b, c, d) representing a^4 + b^4 + c^4 = d^4 
    to sorted list of four (height, mn).
    """
    # Identify odd (c) and even (a, b) terms
    odds = [x for x in abcd[:3] if x % 2 != 0]
    assert 1 == len(odds)
    c = odds[0]
    evens = [x for x in abcd[:3] if x % 2 == 0]
    assert 2 == len(evens)
    a, b= evens[0], evens[1]
    d = QQ(abcd[-1])
    r_raw, s_raw, t = a/d, b/d, c/d

    # Permute (r,s) and flip r, s, t signs.
    h_mn_set = set()
    for r, s, in ((r_raw, s_raw), (s_raw, r_raw),):
        for sr in (r, -r):
            for ss in (s, -s):
                x = (sr + ss) / 2
                y = (sr - ss) / 2
                for st in (t, -t):
                    # Solve Elkies Equation (3c) for u = m/n
                    # Equation: +/- (2u^2 + 1)t^2 = 4(2u^2 - 1)x^2 + 8ux + (1 - 2u^2)
                    # Select + version
                    # Rearranged as Au^2 + Bu + C = 0:
                    # u^2 * [2t^2 - 8x^2 + 2] + u * [-8x] + [t^2 + 4x^2 - 1] = 0
                    quad_A = 2*st**2 - 8*x**2 + 2
                    quad_B = -8*x
                    quad_C = st**2 + 4*x**2 - 1
                    disc = quad_B**2 - 4*quad_A*quad_C
                    if not disc.is_square(): continue
                    root_disc = disc.sqrt()
                    for rd in (root_disc, -root_disc):
                        k = (-quad_B + rd) / (2 * quad_A)
                        # Verify solution with Equation (3b)
                        # (2k^2 + 1)y^2 + (6k^2 - 8k + 3)x^2 + 2(2k^2 - 1)x + 2k = 0
                        check = (2*k**2 + 1)*y**2 + (6*k**2 - 8*k + 3)*x**2 + 2*(2*k**2 - 1)*x + 2*k
                        if 0 != check: continue
                        assert k.numerator()%4 == 0 and k.denominator()%2 == 1
                        h_mn_set.add((k.height(), k))
    return sorted(h_mn_set)

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
    #E_jmin = EllipticCurve(j=E_orig.j_invariant()) # often shorter, but quadratic twist
    E_min = E_orig.minimal_model().short_weierstrass_model()
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
                           coeff_lim: int=3, k_h_lim: int=int(1e15), 
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
    if verbose: print(f'Restrict to {len(ks)} k with height < {k_h_lim:_}')
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

# unfinished
def scan_known_mn(known_s: int=1, known_e: int=9999, max_mn_h: int=100_000_000, 
            quad_h: int=int(1e10), quad_lim: int=5,
            verbose=False, max_d: int=int(1e27)) -> dict:
    """For range of known abcd solutions, generate mn.
    For each mn with height less than max_mn_h:
        Try quad_lim quad_pairs with height less than quad_h.


    ?????
        map u -> d, E a4, a6, E conductor, paired v
            Walk E with v for solutions
    Report by u: #Elliptic Curves, #Conductors
    Return u_into.
    """
    from modular.solutions_modular import get_quartic_pts, \
        u_to_quartic, check_quartic
    for val in list(known.values())[known_s - 1: known_e]:
        abcd = val['abcd']
        if verbose: print(abcd[-1])
        h_mn = abcd_to_h_mn(abcd)
        for h, mn in h_mn:
            assert check_yt(mn)
            if h > max_mn_h: continue
            if verbose: print(f'Trying mn {mn} with height {h}', flush=True)
            quad_hxy = get_rational_points(mn, result_limit=quad_lim)
            for hxy in quad_hxy:
                if hxy[0] > quad_h: break
                quad_xy = hxy[1:]
                q_res = make_quartic(mn, quad_xy)
                if 0 != check_quartic(q_res[2]) or \
                    0 != check_quartic(-q_res[2]): continue
                if verbose:
                    print(f'\tTrying quad_xy {quad_xy}', flush=True)
"""
>>> mn = QQ(-568)/1005
>>> quad_pairs = get_optimized_rational_points(mn)
>>> quad_pairs
[(-106198/581571, 125216/193857), (8225745/54003434, 42633465/54003434), (479811649/1620243202, 909662401/1620243202), (4785497/1806359526, 500466391/602119842), (444634666/1929013011, 447960760/643004337)]
>> y2_coeffs = mn_to_xyt_conics(mn)[0]
>>> y2_coeffs
(1655273, -9532539, 729554, 1141680)
>>> a0, a, b, c = y2_coeffs
>> var('y_var, x_var, k')
(y_var, x_var, k)
>>> eq = a0*y_var**2 - (a*x_var**2 + b*x_var + c)
>>> eq
9532539*x_var^2 + 1655273*y_var^2 - 729554*x_var - 1141680
>>> type(eq)
<class 'sage.symbolic.expression.Expression'>
>>> p = pari(eq).hyperellratpoints(int(1e8), 0)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "cypari2/auto_gen.pxi", line 14712, in cypari2.gen.Gen_base.hyperellratpoints
  File "cypari2/handle_error.pyx", line 211, in cypari2.handle_error._pari_err_handle
cypari2.handle_error.PariError: incorrect type in hyperellratpoints (t_POL)
>>> eq.univariate_polynomial()
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "sage/structure/element.pyx", line 495, in sage.structure.element.Element.__getattr__ (build/cythonized/sage/structure/element.c:13068)
  File "sage/structure/element.pyx", line 508, in sage.structure.element.Element.getattr_from_category (build/cythonized/sage/structure/element.c:13178)
  File "sage/cpython/getattr.pyx", line 363, in sage.cpython.getattr.getattr_from_other_class (build/cythonized/sage/cpython/getattr.c:4677)
AttributeError: 'sage.symbolic.expression.Expression' object has no attribute 'univariate_polynomial'. Did you mean: '_evaluate_polynomial'?

mn = QQ(-568)/1005
a0, a, b, c = y2_coeffs = mn_to_xyt_conics(mn)[0]
R = PolynomialRing(QQ, 'k, y')
k, y = R.gens()
y2 = c/a0 + k * b/a0 + k**2 * a/a0
y2 = y2.univariate_polynomial()
p = pari(y2).hyperellratpoints(int(1e6), 0)
lp = list(p)
pts = [(QQ(k), QQ(y)) for (k, y) in lp]

scan_known_mn(11,111,verbose=True)
...
From here on they take an extremely long time
Investigate get_optimized_rational_points $$$$$$$$$$$$$$$$$$$$$
F26969_14 inx=42, m=(-1873, -200)
26969608212297
Trying mn 200/1873 with height 1873
27497822498977
Trying mn -3696/857 with height 3696
	Trying quad_xy (373281961/4263892434, 135257/781362)
	Trying quad_xy (183552194/7839795627, 165424/390681)
29999857938609
Trying mn -3500/4209 with height 4209
	Trying quad_xy (2524841/130577726, 108429517/130577726)
37352008459537
Trying mn 307700/1565217 with height 1565217
45556888578449
Trying mn -107368/3333 with height 107368
	Trying quad_xy (-547024702/7676412239, 3046936840/7676412239)
"""




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

