"""
Use the uv and their D in solutions_uv.csv with Tomita's elliptic curves
to generates solution.
See my notes in curves_tomita.{md/py}
and Tomita's notes in http://www.maroon.dti.ne.jp/fermat/dioph4e.html
"""
import csv
from multiprocessing import Process, Queue
from sage.all import QQ, ceil, divisors, gcd, lcm, \
        parallel, solve, sqrt, var, \
        Conic, DiagonalQuadraticForm, hilbert_symbol, PolynomialRing

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

def mn_to_xyt_conics(mn: tuple) -> tuple:
    """ Apply rational (m, n) to conics equations in x, y, t
        (2m^2+n^2)y^2 = -(6m^2-8mn+3n^2)x^2 -2(2m^2-n^2)x -2mn        (1) 
        (2m^2+n^2)t^2 = 4(2m^2-n^2)x^2      +8mnx         +(n^2-2m^2) (2)
    Return ((4 coeffs for y^2 eq 1),  (4 coeffs for t^2 eq 2))
    """
    m, n = mn
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
python -m elliptical.solutions_curves mn_to_xyt_conics '(20, -9)'
((881, -4083, -1438, 360), (881, 2876, -1440, -719))
Try other small u associated with Frye solution
>>> mn_to_xyt_conics(((1000, 47)))
((2002209, -5630627, -3995582, -94000), (2002209, 7991164, 376000, -1997791))
>>> mn_to_xyt_conics((-1041, 320))
((2269762, -9474246, -4129924, 666240), (2269762, 8259848, -2664960, -2064962))
>>> mn_to_xyt_conics((-1425, 412))
((4230994, -17389782, -7783012, 1174200), (4230994, 15566024, -4696800, -3891506))
>>> mn_to_xyt_conics((-4209, 3500))
((47681362, -260896086, -46362724, 29463000), (47681362, 92725448, -117852000, -23181362))
>>> mn_to_xyt_conics((30080, 6007))
((1845696849, -4091566067, -3547057502, -361381120), (1845696849, 7094115004, 1445524480, -1773528751))

Try largest v paired with (20, -9)
>>> mn_to_xyt_conics((-6_899_820_729, 369_596_780))
((95351653964462551282, -306456174085512874806, -190156900809779628964, 5100303048031305240), 
 (95351653964462551282, 380313801619559257928, -20401212192125220960, -95078450404889814482))
"""

def is_locally_solvable(mn: tuple) -> bool:
    """ Check real solvability: can ax^2 + bx + c be positive?
    """
    _, a, b, c = y2_coeffs = mn_to_xyt_conics(mn)[0]
    disc = b**2 - 4*a*c
    
    # If a < 0 and disc < 0, RHS is always negative. 
    # Since LHS (y**2 * pos) is always positive, no real solutions.
    return a >= 0 or disc >= 0
"""
>>> is_locally_solvable((20, -9))
True
>>> is_locally_solvable((-9, -20))
False
"""

# Works, but very slow
def find_rational_points_on_conic(mn: tuple, limit: int = 1000):
    # Define Projective Space
    R = PolynomialRing(QQ, 'x, y, z')
    x, y, z = R.gens()
    
    # Coefficients: a0*y^2 = a*x^2 + b*x + c
    a0, a, b, c = mn_to_xyt_conics(mn)[0]
    
    # Homogenize: multiply terms to make every term degree 2
    # a0*y^2 = a*x^2 + b*x*z + c*z^2
    poly = a0*y**2 - (a*x**2 + b*x*z + c*z**2)
    
    C = Conic(poly)
    
    # Local check (p-adic and real)
    if not C.has_rational_point(algorithm='local'):
        print(f'{mn} not locally solvable')
        return None 

    # Find a base point using the 'limit' (bound on coordinate height)
    points = C.rational_points(bound=limit)
    if 0 == len(points):
        print(f'{mn} no rational points found with limit {limit}')
        return None
    # Restrict to positive branch
    return sorted([(x0/z0, y0/z0) for x0, y0, z0 in points if y0 >= 0])

"""
C.rational_points(bound=354)
[(1/354, 75/118), (49/318, 23/106)]

x0, y0, z0 = point # Projective coordinates
    rat_x, rat_y = x0/z0, y0/z0
    
    # 4. Rational Parameterization
    # Once we have one point, we can generate all others using a line of slope 'k'
    t = var('t')
    # Parameterize via the point found
    param = C.rational_parameterization()
    
    return {
        "base_point": (rat_x, rat_y),
        "parameterization": param
    }
"""

def check_solvability(D, delta):
    """ Check the Hilbert symbol on equation X^2 - Dy^2 - delta*z^2 = 0
    at infinite prime (Real Solvability)
    and at 2 and odd primes dividing coeffs
    """
    return hilbert_symbol(D, delta, -1) == 1 == hilbert_symbol(D, delta, 2)

# Works fast
def find_points_by_substitution(mn, limit=1000) -> list:
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
            v = QQ(num) / den
            rhs = D * v**2 + delta
            if rhs >= 0 and rhs.is_square():
                Y = rhs.sqrt()
                #for sign in (1, -1):
                for sign in (1,): # only positive branch
                    x_val = (sign * Y - b) / (2 * a)
                    assert a0 * v *v == a * x_val**2 + b * x_val + c
                    points.append((x_val, v))
    if 0 == len(points):
        print(f'{mn} no rational points found with limit {limit}')
    return points
"""
find_points_by_substitution((20,-9), 107)
[(-73039/144266, -23/106), (49/318, -23/106), (-73039/144266, 23/106), (49/318, 23/106)]
"""

def find_point_instantly(mn: tuple) -> tuple:
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
find_point_instantly((20, -9))
first -1687/5406 -1231/1802
flipped -98423/2452522 -1231/1802
(-1687/5406, -1231/1802)
"""

def get_optimized_rational_points(mn: tuple, 
            search_limit: int=1000, slope_limit: int=50) -> list:
    """ Get many points sorted by height from base, substitution, slopes.
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
        if rhs >= 0 and rhs.is_square():
            y_val = rhs.sqrt()
            # We add both +y and -y
            #for s in (1, -1):
            for s in (1,): # only positive branch
                curr_y = s * y_val
                height = max(x_val.numerator().abs(), x_val.denominator(), 
                             curr_y.numerator().abs(), curr_y.denominator())
                final_points.append((height, x_val, curr_y))

    # Sort by height
    final_points.sort(key=lambda x: x[0])
    return [(p[1], p[2]) for p in final_points]
"""
get_optimized_rational_points((20, -9))[:10]
[(49/318, 23/106), (1/354, 75/118), (-287/598, 211/598), (-63/598, 435/598), 
(94/679, 208/679), (-367/1358, 971/1358), (142/1407, 208/469), 
(126/1447, 696/1447), (185/2058, 325/686), (-1063/3098, 2015/3098)]
"""

def parameterize_conic_solution(mn: tuple, quad_xy: tuple):
    """ Parameterize the conic from a simple point.
    """
    a0, a, b, c = mn_to_xyt_conics(mn)[0]
    x0 = QQ(quad_xy[0][0]) / quad_xy[0][1]
    y0 = QQ(quad_xy[1][0]) / quad_xy[1][1]
    
    # Use Sage's symbolic variables
    var('x_var, k')
    line = k*(x_var - x0) + y0
    
    # The conic equation: a0*y^2 - a2*x^2 - a4*x*z - a6*z^2 = 0
    # For affine x, y (z=1): a0*y^2 - a2*x^2 - a4*x - a6 = 0
    eq = a0*(line**2) - (a*x_var**2 + b*x_var + c)
    
    # Solve for x_var
    roots = solve(eq == 0, x_var)
    
    # Pick the root that isn't the constant x0
    if roots[0].rhs() == x0:
        new_x = roots[1].rhs()
    else:
        new_x = roots[0].rhs()
    
    # Plug new_x into the line equation to get new_y
    new_y = k*(new_x - x0) + y0
    
    # Simplify for a clean output
    return new_x.simplify_full(), new_y.simplify_full()
"""
parameterize_conic_solution((20, -9), ((49, 318), (23, 106)))
(1/318*(43169*k^2 - 121578*k - 657351)/(881*k^2 + 4083), -1/106*(20263*k^2 + 285806*k - 93909)/(881*k^2 + 4083))
"""

# UNTESTED AND LIKELY TO BE SLOW
def find_quartic_points_fast(quartic_poly, den_bound=200):
    # quartic_poly is your f4(k)
    # 1. Define the Hyperelliptic Curve y^2 = f4(k)
    C = HyperellipticCurve(quartic_poly)
    
    # 2. Use the 'ratpoints' algorithm
    # 'bound' is the maximum value for the denominator
    # This is an optimized C-library call
    print(f"Sieving for points with denominator up to {den_bound}...")
    
    # rational_points returns a list of (k, t) pairs
    points = C.rational_points(bound=den_bound, algorithm='ratpoints')
    
    return points

# Example with your coefficients
# R.<k> = QQ[]
# f4 = 4858767860*k^4 - ...
# pts = find_quartic_points_fast(f4, 100)

# TEMPLATE FOR PARALLEL
def persistent_worker(worker_id: int, task_queue, result_queue, inner_loop_gen,
                      chunk_range: tuple, coeffs: tuple):
    # This is the 'Local State' for this specific process
    while True:
        task_range = task_queue.get()
        if task_range is None: # The signal to stop
            result_queue.put(None)
            break
        for r in inner_loop_gen(task_range, coeffs):
            result_queue.put(r)

# TEMPLATE FOR PARALLEL
def chunk_generator(first_v: int, last_v: int, max_chunk: int=50, workers: int=8):
    """ Yield chunks of work to parallel workers.
    Ensures all workers get work by tapering down at the end.
    Allocates from end to beginning because larger could take longer.
    """
    current = last_v

    # Cap chunk size so one chunk doesn't eat the remainder.
    while current >= first_v:
        # Don't let a single chunk take more than remaining work per worker
        # and never exceed a reasonable max (like 50)
        remaining = (current - first_v) + 1
        chunk = max(1, min(remaining // workers, max_chunk))
        yield (current, current - chunk)
        current -= chunk

# TEMPLATE FOR PARALLEL
def find_points_in_y2_conic_gen(chunk_range: tuple, y2_coeffs: tuple):
    """Internal loop for parallel search
    Ignore points on negative branch of isqrt
    """
    a0, a, b, c = y2_coeffs
    bound_e, bound_s = chunk_range
    assert bound_e >= bound_s

    for nz in range(bound_s, bound_e + 1):
        for nx in range(-nz, nz + 1):
            # Using z=nz and x=nx
            # a0*y^2 = a2*nx^2 + a4*nx*nz + a6*nz^2
            rhs = a*nx**2 + b*nx*nz + c*nz**2
            y_sq = QQ(rhs) / a0
            
            if y_sq.is_square():
                ny_val = y_sq.isqrt()
                # We found a point (nx : ny_val : nz)
                yield ((nx, ny_val, nz))

# TEMPLATE FOR PARALLEL
def find_small_point_in_y2_conic_p(mn: tuple, n_cpus: int=8) -> tuple:
    """ Generate y^2 conic for (m, n). Parallel search for rational points.
    Return simplest point.
    """
    y2_coeffs = mn_to_xyt_conics(mn)[0]
    task_queue = Queue()
    result_queue = Queue()

    processes = []
    for i in range(n_cpus):
        p = Process(target=persistent_worker,
                    args=(i, task_queue, result_queue, # Pass both queues
                          y2_coeffs)) # $$$$$$$$$$$ NEED EXPANDING RANGE $$$$$$$$
        p.start()
        processes.append(p)
        
    for chunk in chunk_generator(n_cpus): # $$ NEED EXPANDING RANGE $$$$$$$$
        task_queue.put(chunk)
        
    for _ in range(n_cpus):
        task_queue.put(None)

    # Listen for results while workers are running
    finished_workers = 0
    while finished_workers < n_cpus:
        result = result_queue.get() 
        if result is None:
            finished_workers += 1
        else:
            yield result # Send (u, v, D) back to the main function

    # Cleanup
    for p in processes:
        p.join()


def DEBUG():
    import pdb; pdb.set_trace()
    pass; pass; pass # opportumity to debug

"""
"""

if __name__ == "__main__":
    import sys
    command = eval(sys.argv[1])
    print(command)
    args = list(map(eval, sys.argv[2:]))
    result = command(*args)
    print(result)

""" TODO
Port derive_tomita_quartic
Try find_quartic_points_fast
Rewrite find_rational_point_on_quartic for dynamic parallel.

Look at latest suggestions:
def find_quartic_points_fast(quartic_poly, den_bound=200):
    # quartic_poly is your f4(k)
    # 1. Define the Hyperelliptic Curve y^2 = f4(k)
    C = HyperellipticCurve(quartic_poly)
    
    # 2. Use the 'ratpoints' algorithm
    # 'bound' is the maximum value for the denominator
    # This is an optimized C-library call
    print(f"Sieving for points with denominator up to {den_bound}...")
    
    # rational_points returns a list of (k, t) pairs
    points = C.rational_points(bound=den_bound, algorithm='ratpoints')
    
    return points
# Example with your coefficients
# R.<k> = QQ[]
# f4 = 4858767860*k^4 - ...
# pts = find_quartic_points_fast(f4, 100)



# Convert to Elliptic Curve even without a base point
# Sage uses a more general algorithm (like search for a point first)
E = EllipticCurve_from_hyperelliptic_quartic(f4)

# 'Simplify' the curve (Minimization and Reduction)
E_min = E.minimal_model()
print(f"Minimized Curve: {E_min}")

# Find points on the small curve
# This is much more likely to find 'hidden' points
pts_min = E_min.integral_points(bound=1000)

# You can then map these back to your k-value


# Assuming f4 is your quartic polynomial
from sage.schemes.hyperelliptic_curves.jacobian_generic import HyperellipticJacobian_generic

# This creates the Weierstrass curve associated with your quartic
E = EllipticCurve(f4).jacobian().equivalent_curve()
print(f"Jacobian: {E}")


def fast_sieve_search(f4, limit=200):
    # Pre-calculate f4 mod small primes to eliminate 95% of candidates
    primes = [3, 5, 7, 11, 13, 17, 19, 23]
    squares = {p: set([ (x**2) % p for x in range(p) ]) for p in primes}
    
    for den in range(1, limit):
        for num in range(-limit, limit):
            # The "Sieve Trick": 
            # If (f4(num/den) mod p) isn't a quadratic residue, 
            # it CANNOT be a rational square.
            possible = True
            for p in primes:
                # Use modular arithmetic to keep numbers small
                val_mod_p = (f4.subs(k=QQ(num)/den).numerator()) % p
                # (Need to account for the denominator's inverse mod p too)
                if val_mod_p not in squares[p]:
                    possible = False
                    break
            
            if possible:
                # Only now do we do the expensive square root check
                res = f4.subs(k=QQ(num)/den)
                if res >= 0 and res.is_square():
                    return QQ(num)/den
    return None


    def find_quartic_point_via_minimization(f4):
    # 1. Define the Curve from your Quartic
    C = HyperellipticCurve(f4)
    
    # 2. Map it to an Elliptic Curve (Jacobian)
    # This finds the Weierstrass form and the map 'phi' to get back
    E, phi = C.elliptic_curve(return_morphism=True)
    
    # 3. Minimize the Curve
    # This is the "Magic" step. It reduces your 10^11 coefficients 
    # to the smallest possible integers.
    E_min = E.minimal_model()
    iso = E.isomorphism_to(E_min)
    
    print(f"Original E: {E}")
    print(f"Minimized E: {E_min}")
    
    # 4. Find the Generator
    # On the minimized curve, the points have much smaller 'naive' height.
    # We use Simon's 2-descent which is robust for large curves.
    gens = E_min.simon_two_descent()
    
    if not gens:
        return "No rational points of infinite order found."
    
    # 5. Map the point back to the Quartic
    # Point on E_min -> Point on E -> Point on C (the Quartic)
    P_min = gens[0]
    P_E = iso.inverse()(P_min)
    P_C = phi(P_E)
    
    k_solution = P_C[0]
    return k_solution

# Usage:
# R.<k> = QQ[]
# f4 = 4858767860*k^4 - ...
# print(find_quartic_point_via_minimization(f4))


"""