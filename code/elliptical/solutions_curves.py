"""
Use the u, v and their D in solutions_uv.csv with Tomita's elliptic curves
to generate solutions efficiently.
See log of outputs in solutions_curves.log
See my earlier notes in curves_tomita.{md/py}
and Tomita's notes in http://www.maroon.dti.ne.jp/fermat/dioph4e.html
"""
import csv
from multiprocessing import Process, Queue
from sage.all import Expression, Integer, QQ, Rational, RR, \
    DiagonalQuadraticForm, PolynomialRing, \
    continued_fraction, cos, gcd, hilbert_symbol, lcm, \
    parallel, pi, sage_eval, solve, var
from sage.rings.polynomial.polynomial_rational_flint import Polynomial_rational_flint
from sage.schemes.generic.point import SchemePoint
from time import time
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

def is_locally_solvable(mn: Rational) -> bool:
    """ Check real solvability: can ax^2 + bx + c be positive?
    """
    _, a, b, c = y2_coeffs = mn_to_xyt_conics(mn)[0]
    disc = b**2 - 4*a*c
    
    # If a < 0 and disc < 0, RHS is always negative. 
    # Since LHS (y**2 * pos) is always positive, no real solutions.
    return a >= 0 or disc >= 0
"""
is_locally_solvable(QQ(20/-9))
True
is_locally_solvable(QQ(9/20))
False
"""

def check_solvability(D, delta):
    """ Check the Hilbert symbol on equation X^2 - Dy^2 - delta*z^2 = 0
    at infinite prime (Real Solvability)
    and at 2 and odd primes dividing coeffs
    """
    return hilbert_symbol(D, delta, -1) == 1 == hilbert_symbol(D, delta, 2)

# Works fast
def find_points_by_substitution(mn: Rational, limit: int=1000) -> list:
    a0, a, b, c = mn_to_xyt_conics(mn)[0]
    
    # Y^2 = D*v^2 + delta
    D = 4 * a * a0
    delta = b**2 - 4 * a * c
    if not check_solvability(D, delta):
        print(f'{mn} not locally solvable')
        return None 

    # Search rational v = n/d
    points = []
    for den in range(1, limit):
        for num in range(-limit, limit):
            if gcd(num, den) != 1: continue
            v = QQ(num/den)
            rhs = D * v**2 + delta
            if rhs.is_square():
                Y = rhs.sqrt()
                for sign in (1, -1):
                #for sign in (1,): # only positive branch
                    x_val = (sign * Y - b) / (2 * a)
                    assert a0 * v *v == a * x_val**2 + b * x_val + c
                    points.append((x_val, v))
    if 0 == len(points):
        print(f'{mn} no rational points found with limit {limit}')
    return points
"""
python -m elliptical.solutions_curves find_points_by_substitution -20/9 107
[(-73039/144266, -23/106), (49/318, -23/106), (-73039/144266, 23/106), (49/318, 23/106)]
"""

def find_point_instantly(mn: Rational) -> tuple:
    """ find base point by Gauss-Legendre
    """
    a0, a, b, c = mn_to_xyt_conics(mn)[0]
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
    return (x_val, y_val)
"""
python -m elliptical.solutions_curves find_point_instantly -20/9
(-1687/5406, -1231/1802)
"""

def get_optimized_rational_points(mn: Rational, search_limit: int=1000, 
                    slope_limit: int=50, result_limit: int=10) -> list:
    """ Get small points sorted by height from base, substitution, slopes.
    """
    # Get y2 coefficients
    a0, a, b, c = mn_to_xyt_conics(mn)[0]

    # Find 
    base = find_point_instantly(mn)
    if not base: 
        return []
    bases = find_points_by_substitution(mn, search_limit)
    bases.append(base)
    
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
    print(f'found {len(final_points)}, height {final_points[0][0]} to {final_points[-1][0]}')
    return [(p[1], p[2]) for p in final_points[:result_limit]]
"""
time python -m elliptical.solutions_curves get_optimized_rational_points -20/9
[(49/318, 23/106), (1/354, 75/118), (-287/598, 211/598), (-63/598, 435/598), 
(94/679, 208/679), (-367/1358, 971/1358), (142/1407, 208/469), 
(-391/1418, 1009/1418), (126/1447, 696/1447), (52/1669, 992/1669)]
cpu 47.900 total
"""

def make_quartic(mn: Rational, quad_xy: tuple):
    """ Parameterize y^2 conic with quad point to make quartic
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
    a0_t, a_t, b_t, c_t = t2_coeffs
    t_sq = (a_t*x_k**2 + b_t*x_k + c_t) / a0_t
    
    # 4. Extract the Numerator Polynomial
    # We simplify first to ensure we aren't carrying redundant terms
    poly_expr = t_sq.simplify_full().numerator()
    
    # Convert symbolic expression to a formal polynomial to get clean coefficients
    R = PolynomialRing(QQ, 'k')
    poly = R(poly_expr)
    
    # Handle the GCD and Square Factor Caveat
    coeffs = poly.coefficients()
    common_gcd = gcd([Integer(c) for c in coeffs])
    
    # We divide by common_gcd to get the "Primitive" polynomial
    # Note: If common_gcd is not a square, then Y in Y^2 = clean_poly
    # will be a multiple of the original t by sqrt(common_gcd)
    clean_poly = poly / common_gcd
    
    # 6. Ensure the leading coefficient is positive (Standard Form)
    if clean_poly.leading_coefficient() < 0:
        clean_poly = -clean_poly

    # Return parameterized conic x and y
    # and rhs of y^2 = quartic polynomial, and gcd used to reduce terms
    return x_k, y_k, clean_poly, common_gcd
"""
make_quartic(QQ(20/-9), (QQ(49/318), QQ(23/106)))
(1/318*(43169*k^2 - 121578*k - 657351)/(881*k^2 + 4083), 
 -1/106*(20263*k^2 + 285806*k - 93909)/(881*k^2 + 4083), 
 4858767860*k^4 - 1337905101*k^3 + 32584720500*k^2 - 48737893941*k - 89364400362, 
 4)
"""

def find_quartic_points(quartic_poly: Polynomial_rational_flint, 
            range_s: int=1, range_e: int=1000) -> list:
    """ Search for rational k = nx/nz that make poly(k) or -poly(k) a square.
    primes were selected for those actually used to reject candidates.
    """
    c0, c1, c2, c3, c4 = quartic_poly.list()
    # Search rationals for square rhs
    found = set()
    for nz in range(range_s, range_e + 1):
        nz2 = nz*nz
        nz3 = nz*nz2
        c3nz = c3*nz
        c2nz2 = c2*nz2
        c1nz3 = c1*nz3
        c0nz4 = c0*nz*nz3
        for nx in range(-range_e, range_e + 1):
            if gcd(nx, nz) != 1: continue
            rhs = c4*nx**4 + c3nz*nx**3 + c2nz2*nx**2 + c1nz3*nx + c0nz4
            y2 = abs(rhs)
            if y2.is_square():
                p = QQ(nx/nz)
                print(f'found {p} with rhs {rhs}')
                rhs = quartic_poly(p)
                y2 = abs(rhs)
                if y2.is_square():
                    print(f'\tConfirmed rhs {rhs}')
                    found.add(p)
    return sorted(found)
"""
q_res = make_quartic(QQ(20/-9), (QQ(49/318), QQ(23/106)))
>>> start=time(); find_quartic_points(q_res[2], 1, 1000); elapsed=time()-start; elapsed
found -59/81 with rhs -1493337137920714564
	Confirmed rhs -1493337137920714564/43046721
[-59/81]
5.323254108428955
"""

# slower that find_quartic_points
def find_quartic_point_adaptive(quartic_poly: Polynomial_rational_flint,
            max_denom: int=10_000, n_nodes: int=100, pts_per_node: int=100,
            abs=abs, cos=cos, pi=RR(pi))->list:
    """
    Find rational k that make poly(k) or -poly(k) a square.
    Samples the interval between roots using adaptive density 
    based on the absolute value of the polynomial's derivative.
    """
    coeffs = quartic_poly.list()
    c0, c1, c2, c3, c4 = coeffs
    roots = sorted(quartic_poly.real_roots())
    print(f'roots {roots}')
    assert 2 == len(roots)
    r_min, r_max = roots
    deriv = quartic_poly.derivative()
    
    # Create initial segments as Chebyshev nodes so more dense near roots
    nodes = [0.0] * n_nodes
    for i in range(1, n_nodes + 1):
        # Chebyshev nodes on [-1, 1]
        node = cos((2*i - 1) * pi / (2 * n_nodes))
        # Map from [-1, 1] to [r_min, r_max]
        mapped_node = 0.5 * (r_max + r_min + (r_max - r_min) * node)
        nodes[i-1] = mapped_node
    nodes.extend(roots)
    nodes.sort()
    
    # Calculate local slope intensity
    slopes = [abs(deriv(x)) for x in nodes]
    total_slope = sum(slopes)
    
    # Distribute samples proportional to slope
    seeds = []
    expand = float(pts_per_node * n_nodes / total_slope)
    for i in range(n_nodes):
        # Number of samples in this segment
        n_pts = max(1, int((max(slopes[i], slopes[i+1]) * expand)))
        step = (nodes[i+1] - nodes[i]) / n_pts
        pts = [nodes[i] + step * p for p in range(n_pts)]
        seeds.extend(pts)
    seeds.append(r_max)
      
    # Process seeds to find rational convergents
    candidate_x = set()
    for seed in seeds:
        cf = continued_fraction(seed)
        convergents = cf.convergents()
        for x_frac in convergents:
            if r_min <= x_frac <= r_max:
                if x_frac.denominator() > max_denom: break
                candidate_x.add(x_frac)
    print(f'Testing {len(candidate_x):_} candidates from {len(seeds):_} seeds')
    
    # Check each candidate x for a square |poly(x)|
    found = set()
    for x in candidate_x:
        nx = x.numerator()
        nz = x.denominator()
        rhs = (c4 * nx**4 + 
               c3 * nx**3 * nz + 
               c2 * nx**2 * nz**2 +
               c1 * nx * nz**3 +
               c0 * nz**4)
        y2 = abs(rhs)
        if y2.is_square():
            print(f'found {x}')
            found.add(x)
    return sorted(found)

def k_to_abcd(mn: Rational, quad_xy: tuple, quartic_x: Rational) -> List[Integer]:
    """ Produce solution (A,B,C,D) from (m,n); parameterized x, y; x on quartic.
    """
    q_res = make_quartic(mn, quad_xy)
    x_k, y_k, clean_poly, _ = q_res

    # Evaluate paramterized quartic at specific rational x
    xv = x_k.subs(k=quartic_x)
    yv = y_k.subs(k=quartic_x)
    print(f'x {xv}, y {yv}')
     
    # Calculate t from second conic
    t2_coeffs = mn_to_xyt_conics(mn)[1]
    t_a0, t_a, t_b, t_c = t2_coeffs
    
    # Ensure we stay in QQ
    t_sq = QQ((t_a*xv**2 + t_b*xv + t_c) / t_a0)
    print(f't_sq {t_sq}')
    assert t_sq.is_square()
    tv = t_sq.sqrt()
    print(f't {tv}')

    # Convert to the Elkies/Tomita variables r, s, t
    rv = xv + yv
    sv = xv - yv
    print(f'r {rv}, s {sv}')
    
    # Clear Denominators
    all_fracs = [rv, sv, tv, QQ(1)]
    common_den = lcm([f.denominator() for f in all_fracs])
    
    # Multiply through to get integers
    A = abs(rv * common_den)
    B = abs(sv * common_den)
    C = abs(tv * common_den)
    D = abs(common_den)
    assert A**4 + B**4 + C**4 == D**4
    print(A, B, C, D)

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
More consistency checks
Try negative branch
Try other m,n
"""