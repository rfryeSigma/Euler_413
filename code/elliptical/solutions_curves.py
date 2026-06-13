"""
Use the u, v and their D in solutions_uv.csv with Tomita's elliptic curves
to generate solutions efficiently.
See log of outputs in solutions_curves.log
See my earlier notes in curves_tomita.{md/py}
and Tomita's notes in http://www.maroon.dti.ne.jp/fermat/dioph4e.html
"""
import csv
from datetime import datetime, timedelta
from pdb import set_trace, runcall
from sage.all import GF, Integer, QQ, Rational, RR, \
    DiagonalQuadraticForm, PolynomialRing, \
    continued_fraction, cos, gcd, hilbert_symbol, lcm, \
    oo, pari, pi, sage_eval, solve, var
from sage.rings.polynomial.polynomial_rational_flint import Polynomial_rational_flint
from sage.schemes.generic.point import SchemePoint
from solutions import known
from timeit import repeat
from typing import Union, List, Tuple

def map_D_to_u(file_name: str='solutions_uv.csv') -> dict:
    """
    Map D from solutions_uv.csv to set of u.
    """
    d_map = dict()
    with open(file_name, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file)
        header = next(reader)
        for row in reader:
            u_n, u_d, v_n, v_d, d = [int(row[i]) for i in range(5)]
            d_map[d] = d_map.get(d, set())
            d_map[d].update({(u_n, u_d), (v_n, v_d)})
            # a side check that u obeys mod 4 conditions for m/n
            n, m = (u_n, u_d) if u_n % 2 else (u_d, u_n)
            assert m%4 == 0, f'{m}, {n} fail mod 4'
            #assert abs(n)%4 == 1, f'{m}, {n} fail mod 4'
            #fails 1000, 47
    return d_map
"""
small_d = [d for d in d_map.keys() if d < 1e27]
len(small_d)
85
kd = {x['abcd'][-1] for x in known.values()}
len(kd&small_d)
85
These large missing d probably needed v > 100_000_000
sorted(kd - small_d)
[77107030404994920297, # T77107_20, inx=58, m=(-125, 92),
 101783028910511968041, # P10178_21, inx=59, m=(-9, 20)
 3037467718844497770129, # T30374_22, inx=66, m=(-5, 44),
 9649219915259253551497, # T96492_22, inx=68, m=(-1376, 705)
 18276027741543869996617, # T18276_23, inx=73, m=(-41, 36)
 24504057146788194291849, # T24504_23, inx=74, m=(-5, 44)
 25590155429668179258633, T25590_23, inx=75, m=(-7752, 205)
 29998124444432653523113, T29998_23, inx=76, m=(-1376, 705)
 120175486227071990769561, T12017_24, inx=78, m=(-9, 20)
 1171867103503245199920081, T11718_25, inx=83, m=(-9, 20)
 19874054816411213708481009, B19874_26, inx=89, m=(-5, 8)
 21291952935426564624339201, T21291_26, inx=90, m=(-9, 20)
 96242977191578497031965033, T96242_26, inx=93, m=(-41, 36)
 133140691304639620846181457, B13314_27, inx=95, m=(-5, 8)
 227529118288906398066378489, P22752_27, inx=97, m=(-9, 20)
]
Only 9 u for 422_481. Other 3 are in missing d
(Pdb) d_map[422481]
{(3521543, 9580960), 
(-9, 20), 
(-1041, 320), 
(52463, 660460), 
(30080, 6007), 
(330353, 48940), 
(1000, 47), 
(-4209, 3500), 
(29957400, 6538471)}

Yes, the missing v for 422_481 all have denom > 100_000_000
du = abcd_to_u_set((95_800, 217_519, 414_560, 422_481))
du = {(numerator(u), denominator(u)) for u in du}
du - d_map[422481]
{(-71490240, 101_943_281), (-167767337, 43_538_900), (-6_899_820_729, 369_596_780)}
"""

def mn_to_xyt_conics(mn: Rational) -> tuple:
    """ Apply m/n to conics equations in x, y, t
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
""""
python -m elliptical.solutions_curves mn_to_xyt_conics 20/-9
((881, -4083, -1438, 360), (881, 2876, -1440, -719))
Try other small u associated with Frye solution
mn_to_xyt_conics(QQ(1000/47))
((2002209, -5630627, -3995582, -94000), (2002209, 7991164, 376000, -1997791))
mn_to_xyt_conics(QQ(-1041/320))
((2269762, -9474246, -4129924, 666240), (2269762, 8259848, -2664960, -2064962))
"""

# unused
def is_locally_solvable(mn: Rational) -> bool:
    """ Check real solvability: can ax^2 + bx + c be positive?
    """
    _, a, b, c = y2_coeffs = mn_to_xyt_conics(mn)[0]
    disc = b**2 - 4*a*c
    # If a < 0 and disc < 0, RHS is always negative. 
    return a >= 0 or disc >= 0
"""
is_locally_solvable(QQ(20/-9))
True
is_locally_solvable(QQ(9/20))
False
"""

def check_quadratic(D, delta):
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
        if not check_quadratic(D, delta):
            return False
        Q = DiagonalQuadraticForm(QQ, [1, -D, -delta])
        try: # find a rational point (X, V, Z)
            point = Q.solve()
        except Exception as e: return False
    return True

def find_point_instantly(mn: Rational) -> tuple:
    """ Find base point by Gauss-Legendre. Check both t^2 and y^2,
    but return point for y^2.
    """
    y2_coeffs, t2_coeffs = mn_to_xyt_conics(mn)
    for a0, a, b, c in (t2_coeffs, y2_coeffs):
        D = 4 * a * a0
        delta = b**2 - 4 * a * c
        
        # We want to solve X^2 - D*V^2 - delta*Z^2 = 0
        # Create the quadratic form matrix:
        # [ 1  0      0   ]
        # [ 0 -D      0   ]
        # [ 0  0 -delta   ]
        Q = DiagonalQuadraticForm(QQ, [1, -D, -delta])
        
        # This finds a rational point (X, V, Z) instantly
        try:
            point = Q.solve()
        except Exception as e:
            print(f"Error: {e}")
            return None
        assert point is not None
            
        X, V, Z = point
        # Convert back to original x, y
        x_val = (X / Z - b) / (2 * a)
        y_val = V / Z
        assert a0 * y_val**2 == a * x_val**2 + b * x_val + c
        # print(x_val, y_val)
    return (x_val, y_val)
"""
python -m elliptical.solutions_curves find_point_instantly -20/9
(-1687/5406, -1231/1802)

find_point_instantly(QQ(20/9))
103780/166089 281/231
Error: no solution found (local obstruction at 7)
find_point_instantly(QQ(9/20))
Error: no solution found (local obstruction at 101)
find_point_instantly(QQ(20/9))
103780/166089 281/231
Error: no solution found (local obstruction at 7)
find_point_instantly(QQ(-9/20))
Error: no solution found (local obstruction at 101)
find_point_instantly(QQ(20/-9))
260/231 281/231
-1687/5406 -1231/1802
(-1687/5406, -1231/1802)
so only 20/-9 has no local obstruction to both t^2 and y^2.
"""

# obsolete
def find_points_by_substitution(mn: Rational, limit: int=200, 
                                first_d: int=1) -> list:
    a0, a, b, c = mn_to_xyt_conics(mn)[0]
    
    # Y^2 = D*v^2 + delta
    D = 4 * a * a0
    delta = b**2 - 4 * a * c
    if not check_quadratic(D, delta):
        print(f'{mn} not locally solvable')
        return None 

    # Search rational v = n/d
    points = []
    for den in range(first_d, limit):
        for num in range(-limit, limit):
            if gcd(num, den) != 1: continue
            v = QQ(num/den)
            rhs = D * v**2 + delta
            if rhs.is_square():
                Y = rhs.sqrt()
                for sign in (1, -1):
                #for sign in (1,): # only positive branch
                    x_val = (sign * Y - b) / (2 * a)
                    assert a0 * v * v == a * x_val**2 + b * x_val + c
                    points.append((x_val, v))
    #if 0 == len(points):
        #print(f'{mn} no quadratic points with limit {limit}')
    return points
"""
python -m elliptical.solutions_curves find_points_by_substitution -20/9 107
[(-73039/144266, -23/106), (49/318, -23/106), (-73039/144266, 23/106), (49/318, 23/106)]
"""

# obsolete
def get_optimized_rational_points(mn: Rational, batch_size: int=200, 
                    slope_limit: int=20, result_limit: int=5) -> list:
    """ Get small points sorted by height from base, substitution, slopes.
    """
    # Get y2 coefficients
    a0, a, b, c = mn_to_xyt_conics(mn)[0]

    # Collect base points by Gauss-Legendre and substitution batches
    base = find_point_instantly(mn)
    if not base: return []
    bases = [base]
    first_d = 1
    for i in range(1, 4):
        batch = find_points_by_substitution(mn, i * batch_size, first_d)
        bases.extend(batch)
        if len(bases) > max(15, 3 * result_limit): break
        first_d += i * batch_size

    # Augment the bases with x-axis reflection and collect for search
    seeds = set()
    for x0, y0 in bases:
        seeds.update([(x0, y0), ((-b - 2*a*x0)/(a), y0)])
 
    # Collect unique x-coordinates via slopes
    unique_x = {s[0] for s in seeds}
    for sx, sy in seeds:
        for den in range(1, slope_limit):
            for num in range(-slope_limit, slope_limit):
                k = QQ(num) / den
                denom = a0 * k**2 - a
                if denom == 0: continue
                
                x1 = (sx * (a + a0 * k**2) + 2 * a0 * k * sy + b) / denom
                unique_x.add(x1)

    # Re-solve for y for unique x to find the 'cleanest' y.
    final_points = []
    for x_val in unique_x:
        rhs = (a * x_val**2 + b * x_val + c) / a0
        if rhs.is_square():
            y_val = rhs.sqrt()
            # We add both +y and -y
            #for s in (1, -1):
            for s in (1,): # only positive branch
                curr_y = s * y_val
                height = max(x_val.height(), curr_y.height())
                final_points.append((height, x_val, curr_y))

    # Sort by height
    final_points.sort(key=lambda x: x[0])
    #print(f'found {len(final_points)}, height {final_points[0][0]} to {final_points[-1][0]}')
    return [(p[1], p[2]) for p in final_points[:result_limit]]
"""
time python -m elliptical.solutions_curves get_optimized_rational_points -20/9
[(49/318, 23/106), (1/354, 75/118), (-287/598, 211/598), (-63/598, 435/598), 
(94/679, 208/679), (-367/1358, 971/1358), (142/1407, 208/469), 
(-391/1418, 1009/1418), (126/1447, 696/1447), (52/1669, 992/1669)]
cpu 47.900 total
"""

def get_rational_points(mn: Rational, max_pt: int=1_000_000,
        result_limit: int=5, verbose=True,
        d_list: tuple=(100_000, 1_000_000, 10_000_000,),
        n_mult: tuple= (5, 2, 1,)) -> list:
    """Get result_limit rational points on y2 polynomial
    """
    # Build y2 polynomial
    a0, a, b, c = mn_to_xyt_conics(mn)[0]
    R = PolynomialRing(QQ, 'k, y')
    k, y = R.gens()
    y2 = c/a0 + k * b/a0 + k**2 * a/a0
    y2 = y2.univariate_polynomial()

    # Collect base points by Gauss-Legendre and hyperellratpoints batches
    base = find_point_instantly(mn)
    if not base: return []
    bases = [base]
    d = 0
    for r_inx in range(len(d_list)):
        s = d + 1
        if max_pt <= d: break
        d = min(d_list[r_inx], max_pt)
        n = n_mult[r_inx] * d
        p = pari(y2).hyperellratpoints([n, [s, d]], 0)
        lp = list(p)
        batch = [(QQ(k), QQ(y)) for (k, y) in lp if 0 <= y]
        bases.extend(batch)
        if len(bases) > 1 * result_limit: break

    # Augment the bases with x-axis reflection and collect for search
    seeds = set()
    for x0, y0 in set(bases):
        seeds.update([(x0, y0), ((-b - 2*a*x0)/(a), y0)])

    # Sort by height
    final_points = []
    for x0, y0 in seeds:
        height = max(x0.height(), y0.height())
        final_points.append((height, x0, y0))
    return sorted(final_points)[:result_limit]
    
def make_quartic(mn: Rational, quad_xy: tuple):
    """ Parameterize y^2 and t^2 conic with quad point to make t^2 quartic
    """
    y2_coeffs, t2_coeffs = mn_to_xyt_conics(mn)
    x0, y0 = quad_xy

    # Use Sage's symbolic variables
    var('x_var, k')
    line = k*(x_var - x0) + y0
    
    # Parameterize affine conic equation
    a0, a, b, c = y2_coeffs
    eq = a0*(line**2) - (a*x_var**2 + b*x_var + c)
    
    # Solve for x_var
    roots = solve(eq == 0, x_var)
    
    # Pick the root that isn't the constant x0
    if roots[0].rhs() == x0:
        new_x = roots[1].rhs()
    else:
        new_x = roots[0].rhs()
    x_k = new_x.simplify_full()
    #print(f'x_k {x_k}')

    # Plug new_x into the line equation to get new_y
    y_k = k*(x_k - x0) + y0
    y_k = y_k.simplify_full()
    #print(f'y_k {y_k}')

    # Use parameterized x_k to calculate affine t^2 symbolically
    t_a0, t_a, t_b, t_c = t2_coeffs
    t_sq = (t_a*x_k**2 + t_b*x_k + t_c) / t_a0
    
    # Extract the Numerator Polynomial and simplify
    poly_expr = t_sq.simplify_full().numerator()
    #print(poly_expr)
    
    # Convert symbolic expression to a formal polynomial to get clean coefficients
    R = PolynomialRing(QQ, 'k')
    try:
        poly = R(poly_expr)
    except Exception as e:
        print(f'Error converting to polynomial: {e}')
        print(f'\tfor mn {mn}, quad_xy {quad_xy}, poly_expr: {poly_expr}')
        return None
    
    # Handle the GCD and Square Factor Caveat
    coeffs = poly.coefficients()
    common_gcd = gcd([Integer(c) for c in coeffs])
    assert common_gcd.is_square()
    clean_poly = poly / common_gcd
    
    # Ensure the leading coefficient is positive (Standard Form)
    if clean_poly.leading_coefficient() < 0:
        clean_poly = -clean_poly

    # Force t_sq into the Fraction Field of the Polynomial Ring
    # This strips away the fragile symbolic 'var' types entirely.
    Frac_R = R.fraction_field()
    t_sq_rational = Frac_R(t_sq) 
    
    # Now divide formal rational expressions and let sage simplify.
    t_scale_expr = t_sq_rational / clean_poly

    # Return parameterized conics x and y
    # and rhs of y^2 = quartic polynomial
    # and how to evaluate t_sq
    return x_k, y_k, clean_poly, t_scale_expr
"""
make_quartic(QQ(20/-9), (QQ(49/318), QQ(23/106)))
(1/318*(43169*k^2 - 121578*k - 657351)/(881*k^2 + 4083), 
 -1/106*(20263*k^2 + 285806*k - 93909)/(881*k^2 + 4083), 
 4858767860*k^4 - 1337905101*k^3 + 32584720500*k^2 - 48737893941*k - 89364400362, 
-4/19622126241/(k^4 + 8166/881*k^2 + 16670889/776161))
"""

# Old. See more robust solutions_modular.get_quartic_pts
def find_quartic_points_hyper(quartic_poly: Polynomial_rational_flint,
            max_num: int, min_den: int, max_den: int) -> list:
    """ Use hyperellratpoints on quartic_poly to find quartic points.
    """
    poly = quartic_poly
    roots = poly.real_roots()
    print(f'quartic roots {roots}')
    assert 2 == len(roots)
    #if 0 < poly.list()[-1]: poly = -poly
    pts = pari(poly).hyperellratpoints((max_num, (min_den, max_den)), 0)
    if 0 < len(pts): return pts
    pts = pari(-poly).hyperellratpoints((max_num, (min_den, max_den)), 0)
    return pts

def k_to_abcd(mn: Rational, quad_xy: tuple, quartic_x: Rational) -> List[Integer]:
    """ Produce solution (A,B,C,D) from (m,n); parameterized x, y; x on quartic.
    """
    q_res = make_quartic(mn, quad_xy)
    x_k, y_k, quartic_poly, t_scale_expr = q_res

    # Check that quartic_poly is square at x
    rhs = quartic_poly.subs(k=quartic_x)
    y2 = abs(rhs)
    assert y2.is_square()

    # Evaluate parameterized conic at x
    xv = x_k.subs(k=quartic_x)
    yv = y_k.subs(k=quartic_x)
     
    # Evaluate the scale factor at our target k
    scale_val = t_scale_expr.subs(k=quartic_x)
    
    # t_sq is the clean_poly value times the scale factor
    t_sq = abs(QQ(y2 * scale_val))
    #print(f't_sq {t_sq}')
    assert t_sq.is_square(), 't_sq has a quadratic twist'
    tv = t_sq.sqrt()
    
    # Convert to the Elkies/Tomita variables r, s, t
    rv = xv + yv
    sv = xv - yv
    #print(f'r {rv}, s {sv}, t {tv}')
    
    ## Clear Denominators
    all_fracs = [rv, sv, tv, QQ(1)]
    common_den = lcm([f.denominator() for f in all_fracs])
    
    # Multiply through to get integers
    A = abs(rv * common_den)
    B = abs(sv * common_den)
    C = abs(tv * common_den)
    D = abs(common_den)
    
    # Verification
    assert A**4 + B**4 + C**4 == D**4

    # Return sorted A,B,C and the D
    lhs_list = sorted([Integer(A), Integer(B), Integer(C)])
    return lhs_list, Integer(D)
"""
k_to_abcd(QQ(20/-9), (QQ(49/318), QQ(23/106)), QQ(-59/81))
x -159380/422481, y 85060/140827
t_sq 47314515361/178490195361
t 217519/422481
r 95800/422481, s -414560/422481
95800 414560 217519 422481
([95800, 217519, 414560], 422481)

mn = QQ(-291)/80
quad_xy = (QQ(-1448)/7293, QQ(152)/221)
k0 = QQ(1292115)/175762
k_to_abcd(mn, quad_xy, k0)
    assert t_sq.is_square(), 't_sq has a quadratic twist'
AssertionError: t_sq has a quadratic twist
t_sq
1053976225/692082399137
The numerator is square, but not the denominator.
I don't understand why this fails.
Note that -40/291 fails with same quad_xy:
>>> k_to_abcd(QQ(-40)/291, (QQ(1448)/7293, QQ(152)/221), QQ(-1292115)/175762)
AssertionError: t_sq has a quadratic twist
"""

# Old. See find_mn_y2t2_pts
def search_mn(m: int, n_s: int, n_e: int, quad_n: int=5, 
              quart_s: int=1, quart_e: int=50_000) -> set:
    """ Given m, search n over range (n_s, n_e)
    for quartics with rational point.
    For each n, try quad_n pts, and for quart, try up to quart_e.
    """
    assert m%4 == 0
    assert n_s%4 == 1
    assert 1 <= n_s <= n_e
    d_dict = dict()
    kd = {v['abcd'][-1] for v in known.values()}
    for n in range(n_s, n_e + 1, 4):
        if gcd(m , n) != 1: continue
        print(f'Trying n {n}')
        for mn_pair in ((m,n), (n,m), (m,-n), (-n,m)):
            mn = QQ(mn_pair[0]) / QQ(mn_pair[1])
            if find_point_instantly(mn) == None: continue
            print(f'Trying mn {mn}')
            quad_hxy = get_rational_points(mn, result_limit=quad_n)
            for hxy in quad_hxy:
                quad_xy = hxy[1:]
                q_res = make_quartic(mn, quad_xy)
                print(f'Trying mn {mn}, quad_xy {quad_xy}')
                assert 2 == len(q_res[2].real_roots())
                pts = find_quartic_points_hyper(q_res[2], quart_s, quart_e)
                if 0 == len(pts):
                    print(f'No quartic point for mn {mn}')
                for pt in pts:
                    _, d = k_to_abcd(mn, quad_xy, pt)
                    if d > 1e27: continue
                    d_dict[d] = d_dict.get(d, 0) + 1
                    assert d in kd
    return d_dict

def collect_mn(m: int, n_s: int, n_e: int, mod: int=4) -> list:
    """ Given m, search n over range (n_s, n_e)
    for quadratics with rational points.
    keep each mn with no local obstruction.
    """
    assert m%mod == 0
    assert n_s%mod == 1
    assert 1 <= n_s <= n_e
    found_mn = []
    for n in range(n_s, n_e + 1, mod):
        if gcd(m , n) != 1: continue
        #print(f'Trying n {n}')
        for mn_pair in ((m,n), (n,m), (m,-n), (-n,m)):
            mn = QQ(mn_pair[0]) / QQ(mn_pair[1])
            if find_point_instantly(mn) == None: continue
            found_mn.append(mn)
    return found_mn
"""
mn_list = collect_mn(20, 1, 10_000, 10)
len(mn_list)
135
collected_mn_to_brute(mn_list)
Counts: mn 135, D 0, big 0, known 0
doesn't find anything
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
