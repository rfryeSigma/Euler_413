"""
Code to explore how to manipulate elliptic curves like Seiji Tomita does.
See my notes in curves_tomita.md
and Tomita's notes in http://www.maroon.dti.ne.jp/fermat/dioph4e.html
"""
from sage.all import Integer, QQ, lcm, ceil, gcd, solve, sqrt, var, \
        EllipticCurve, parallel, PolynomialRing

def mn_to_coeffs_y2(m: int, n: int) -> tuple:
    """Apply m, n to Tomita's eq 1 (Elkies' eq 3b)
    (2m^2+n^2)y^2=-(6m^2-8mn+3n^2)x^2-2(2m^2-n^2)x-2mn
    and return the evaluated coefficients (a0, a2, a4, a6)
    for Tomita's eq. 3
    """
    m2 = m*m
    n2 = n*n
    mn = m*n
    a0 = 2 * m2 + n2
    a2 = -(6 * m2 - 8 * mn + 3 *n2)
    a4 = -2 * (2 * m2 - n2)
    a6 = -2 * mn
    return (a0, a2, a4, a6)
""""
Tomita's eq. 3
% python -m elliptical.curves_tomita mn_to_coeffs_y2 20 -9
(881, -4083, -1438, 360)
"""

def mn_to_coeffs_t2(m: int, n: int) -> tuple:
    """Apply m, n to Tomita's eq 2 (Elkie's eq 3c)
    (2m2+n2)t2=4(2m2-n2)x2+8mnx+(n2-2m2)
    and return the evaluated coefficients (a0, a2, a4, a6)
    for Tomita's eq. 4
    """
    m2 = m*m
    n2 = n*n
    mn = m*n
    a0 = 2 * m2 + n2
    a2 = 4 * (2 * m2 - n2)
    a4 = 8 * mn
    a6 = n2 - 2 * m2
    return (a0, a2, a4, a6)
"""
Tomita's eq. 4
% python -m elliptical.curves_tomita mn_to_coeffs_t2 20 -9
(881, 2876, -1440, -719)
"""

def find_rational_points_in_conic(range_s: int, range_e: int, m: int, n: int,
                                  ncpus: int=8) -> list:
    """ Find rational points on Tomita's eq 3
    """
    @parallel(ncpus=ncpus)
    def find_pts_worker(bound_s, bound_e, coeffs):
        a0, a2, a4, a6 = coeffs
        found_in_worker = []
        
        # Projective equation: a0*y^2 = a2*x^2 + a4*x*z + a6*z^2
        # For x/z, we check numerator 'nx' and denominator 'nz'
        
        # To find points like (1/354), we search denominators up to bound_e
        # Sage's C.rational_points(bound=400) searches HEIGHT.
        # Height of (a/b : c/d : 1) is max(|a*d|, |b*c|, |b*d|)
        
        for nz in range(bound_s, bound_e + 1):
            for nx in range(-nz, nz + 1):
                # Using z=nz and x=nx
                # a0*y^2 = a2*nx^2 + a4*nx*nz + a6*nz^2
                rhs = a2*nx**2 + a4*nx*nz + a6*nz**2
                y_sq = QQ(rhs) / a0
                
                if y_sq >= 0 and y_sq.is_square():
                    ny_val = y_sq.sqrt()
                    # We found a point (nx : ny_val : nz)
                    found_in_worker.append((nx, ny_val, nz))
                    if ny_val != 0:
                        found_in_worker.append((nx, -ny_val, nz))
        return found_in_worker

    a_coeffs = mn_to_coeffs_y2(m, n) # (881, -4083, -1438, 360)
    
    # We split by denominator 'nz' across CPUs
    r = int(ceil((range_e - range_s + 1) / ncpus))
    tasks = [(range_s + i*r, range_s + (i+1)*r - 1, a_coeffs) for i in range(ncpus)]
    
    final_pts = set()
    for task_input, found in find_pts_worker(tasks):
        if isinstance(found, list):
            for nx, ny, nz in found:
                if nz != 0:
                    # Construct coordinates as simplified fractions
                    # This turns (49, 69, 318) into (49/318, 23/106)
                    final_pts.add((QQ(nx)/nz, QQ(ny)/nz))
                else:
                    # Point at infinity (1 : y/x : 0)
                    final_pts.add((QQ(nx), QQ(ny), 0))
    
    return sorted(list(final_pts))
"""
Find rational points on Tomita's eq 3
python -m elliptical.curves_tomita find_rational_points_in_conic 1 1000 20 -9
<function find_rational_points_in_conic at 0x12e622840>
[(-287/598, -211/598), (-287/598, 211/598), (-63/598, -435/598), (-63/598, 435/598), 
(1/354, -75/118), (1/354, 75/118), (94/679, -208/679), (94/679, 208/679), 
(49/318, -23/106), (49/318, 23/106)]
"""

def parameterize_conic_solution(m: int, n: int, x0_num: int, x0_den: int,
                                y0_num: int, y0_den: int):
    a0, a2, a4, a6 = mn_to_coeffs_y2(m, n)
    """ Use a simple rational point on Tomita's eq 3 to parameterize it
    as Tomita's eq. 5 and 6
    """
    x0 = QQ(x0_num) / x0_den
    y0 = QQ(y0_num) / y0_den
    
    # Use Sage's symbolic variables
    var('x_var, k')
    line = k*(x_var - x0) + y0
    
    # The conic equation: a0*y^2 - a2*x^2 - a4*x*z - a6*z^2 = 0
    # For affine x, y (z=1): a0*y^2 - a2*x^2 - a4*x - a6 = 0
    eq = a0*(line**2) - a2*x_var**2 - a4*x_var - a6
    
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
This parameterization agrees with Tomita equation 5, 6.
% python -m elliptical.curves_tomita parameterize_conic_solution 20 -9 49 318 23 106 
(1/318*(43169*k^2 - 121578*k - 657351)/(881*k^2 + 4083), -1/106*(20263*k^2 + 285806*k - 93909)/(881*k^2 + 4083))

These give different, but presumably equivalent results: 
% python -m elliptical.curves_tomita parameterize_conic_solution 20 -9 49 318 -23 106
(1/318*(43169*k^2 + 121578*k - 657351)/(881*k^2 + 4083), 1/106*(20263*k^2 - 285806*k - 93909)/(881*k^2 + 4083))
% python -m elliptical.curves_tomita parameterize_conic_solution 20 -9 1 354 75 118  
(1/354*(881*k^2 - 396450*k - 513135)/(881*k^2 + 4083), -1/118*(66075*k^2 + 172406*k - 306225)/(881*k^2 + 4083))
"""

def derive_tomita_quartic(m, n, x0_num, x0_den, y0_num, y0_den):
    """ Substitute Tomita eq 5 into eq 4 to get parameterized eq 7.
    """
    # 1. Get coefficients for Tomita (4)
    a0, a2, a4, a6 = mn_to_coeffs_t2(m, n)
    
    # 2. Get the x(k) parametrization
    x_k, _ = parameterize_conic_solution(m, n, x0_num, x0_den, y0_num, y0_den)
    
    # 3. Calculate t^2 symbolically
    # a0 * t^2 = a2*x^2 + a4*x + a6
    t_sq = (a2*x_k**2 + a4*x_k + a6) / a0
    
    # 4. Extract the Numerator Polynomial
    # We simplify first to ensure we aren't carrying redundant terms
    k = var('k')
    poly_expr = t_sq.simplify_full().numerator()
    
    # Convert symbolic expression to a formal polynomial to get clean coefficients
    R = PolynomialRing(QQ, 'k')
    poly = R(poly_expr)
    
    # 5. Handle the GCD and Square Factor Caveat
    coeffs = poly.coefficients()
    common_gcd = gcd([Integer(c) for c in coeffs])
    
    # We divide by the common_gcd to get the "Primitive" polynomial
    # Note: If common_gcd is not a square, the Y in Y^2 = clean_poly
    # will be a multiple of the original t by sqrt(common_gcd)
    clean_poly = poly / common_gcd
    
    # 6. Ensure the leading coefficient is positive (Standard Form)
    if clean_poly.leading_coefficient() < 0:
        clean_poly = -clean_poly
        
    return clean_poly, common_gcd
"""
% python -m elliptical.curves_tomita derive_tomita_quartic 20 -9 49 318 23 106
(4858767860*k^4 - 1337905101*k^3 + 32584720500*k^2 - 48737893941*k - 89364400362, 4)

If I multiply my answer by the gcd 4, I get Tomita's eq 7
Y2=19435071440k4-5351620404k3+130338882000k2-194951575764k-357457601448

But in http://www.maroon.dti.ne.jp/fermat/grouplaw1e.html, he gets
V2=-19435071440U4+5351620404U3-130338882000U2+194951575764U+357457601448.......(3)
which he transforms to 
Minimal Weierstrass form: y2 = x3+2265722465761x-3154189403034549278.............(4)
with points
p1=[978559, 0]
p2=[47971729/49, 16603172706/343]
p3=[1237921, 1244044242]
from which he uses the Group Law to solve
    (n1,n2,n3)=(0,1,0):MacLeod's case
    (n1,n2,n3)=(1,0,0):Frye's case
and 4 others
"""

def find_rational_point_on_quartic(m, n, x0_num, x0_den, y0_num, y0_den, range_s: int, range_e: int, ncpus: int=8):
    """ Search for rational k = nx/nz that make poly(k) or -poly(k) a square.
    """
    poly, _ = derive_tomita_quartic(m, n, x0_num, x0_den, y0_num, y0_den)
    # Standard form: c4*k^4 + c3*k^3 + c2*k^2 + c1*k + c0
    print(poly)
    coeffs = poly.list() 
    print(coeffs)

    @parallel(ncpus=ncpus)
    def find_k_worker(nz_start, nz_end, poly_coeffs):
        found_in_worker = []
        c0, c1, c2, c3, c4 = poly_coeffs
        
        for nz in range(nz_start, nz_end + 1):
            nz2 = nz*nz
            nz3 = nz2*nz
            nz4 = nz3*nz
            for nx in range(-nz, nz + 1):
                # We evaluate poly(nx/nz) * nz^4 to keep it in integers:
                # Y_int^2 = c4*nx^4 + c3*nx^3*nz + c2*nx^2*nz^2 + c1*nx*nz^3 + c0*nz^4
                rhs = (c4 * nx**4 + 
                       c3 * nx**3 * nz + 
                       c2 * nx**2 * nz**2 + 
                       c1 * nx * nz**3 + 
                       c0 * nz**4)
                
                # Check both positive and negative (handling the twist)
                if rhs.is_square():
                    found_in_worker.append(QQ(nx)/nz)
                elif (-rhs).is_square():
                    found_in_worker.append(QQ(nx)/nz)
                    
        return found_in_worker

    # Split the search range (nz) across cores
    r = int(ceil((range_e - range_s + 1) / ncpus))
    tasks = [(range_s + i*r, range_s + (i+1)*r - 1, coeffs) for i in range(ncpus)]
    
    found_k_values = set()
    for task_input, results in find_k_worker(tasks):
        for k_val in results:
            found_k_values.add(k_val)
            
    return sorted(list(found_k_values))
"""
The workers are unbalanced; it would be better to use a dynamic scheduler.  $$$$$$$$$$$$$$

% python -m elliptical.curves_tomita find_rational_point_on_quartic 20 -9 49 318 23 106 1 1000
4858767860*k^4 - 1337905101*k^3 + 32584720500*k^2 - 48737893941*k - 89364400362
[-89364400362, -48737893941, 32584720500, -1337905101, 4858767860]
[-59/81]
Over various ranges up to 100_000, only this point shows up.
"""

def verify_step_7(poly_coeffs: list=[4858767860, -1337905101, 32584720500, -48737893941, -89364400362],
                  common_gcd: int=4, k_num: int=-59, k_den: int=81):
    """ Verify that the given k value yields a Y^2.
    """
    k = QQ(k_num) / k_den
    print([common_gcd * c for c in poly_coeffs])
     
    # Calculate RHS: c4*k^4 + c3*k^3 + c2*k^2 + c1*k + c0
    rhs = sum(c * k**i for i, c in enumerate(reversed(poly_coeffs)))
        
    print(f"RHS value: {rhs}")
    if rhs < 0:
        print('RHS is negative so switching sign')
        rhs = - rhs
    print(f"Is RHS a square? {rhs.is_square()}")
    if rhs.is_square():
        print(f"Y = {rhs.sqrt()}")
        return
"""
% python -m elliptical.curves_tomita verify_step_7
[19435071440, -5351620404, 130338882000, -194951575764, -357457601448]
RHS value: -1493337137920714564/43046721
RHS is negative so switching sign
Is RHS a square? True
Y = 1222021742/6561
"""

def model_quartic_as_elliptic_curve(m, n, x0_num, x0_den, y0_num, y0_den, pt_num, pt_den):
    """Convert parameterized quartic to Elliptic Curve
    """
    poly, _ = derive_tomita_quartic(m, n, x0_num, x0_den, y0_num, y0_den)
    print('poly', poly)

    k0 = QQ(pt_num) / pt_den
    print('k0', k0)

    # Calculate the Y coordinate for k0. Twist if negative.
    y_sq = poly.subs(k=k0)    
    if not y_sq.is_square():
        poly = -poly
        y_sq = poly.subs(k=k0)        
    y0 = y_sq.sqrt()
    print('y0', y0, 'y0^2', y0*y0)
   
     # u = k - k0 is point shifted to origin
    R = PolynomialRing(QQ, 'u')
    u = R.gen()
    shifted_poly = poly.subs({poly.parent().gen(): u + k0})
    coeffs = shifted_poly.list()
    
    # Standard form: Y^2 = a*u^4 + b*u^3 + c*u^2 + d*u + e^2
    # Ensure we have 5 coefficients
    while len(coeffs) < 5: coeffs.append(0)
    e2, d, c, b, a = coeffs
    print('Elliptic coeffs with e2\n', coeffs)
    e = y0 # The sqrt of the constant term
    
    # 2. Coefficients of the Weierstrass form y^2 + a1*xy + a3*y = x^3 + a2*x^2 + a4*x + a6
    # These formulas come from the Mordell/Cassels transformation
    a1 = d/e
    a2 = c - (d**2)/(4*e**2)
    a3 = 2*e*b
    a4 = -4*e**2 * a
    a6 = a2 * a4
    
    E = EllipticCurve([a1, a2, a3, a4, a6])
    E_short = E.short_weierstrass_model()
    E_min = E_short.minimal_model()
    E_jmin = EllipticCurve(j=E.j_invariant()) # often shorter, but quadratic twist
    return E_min, E_jmin

"""
% python -m elliptical.curves_tomita model_quartic_as_elliptic_curve 20 -9 49 318 23 106 -59 81
poly 4858767860*k^4 - 1337905101*k^3 + 32584720500*k^2 - 48737893941*k - 89364400362
k0 -59/81
y0 1222021742/6561 y0^2 1493337137920714564/43046721
Elliptic coeffs with e2
 [1493337137920714564/43046721, 56251607253382924/531441, -111483374052499/2187, 1255039528141/81, -4858767860]
(Elliptic Curve defined by y^2 = x^3 + 2265722465761*x - 3154189403034549278 over Rational Field, 
Elliptic Curve defined by y^2 = x^3 + 2265722465761*x - 3154189403034549278 over Rational Field)

which matches C1:X^3+2265722465761X-3154189403034549278=Y^
Check points used in http://www.maroon.dti.ne.jp/fermat/grouplaw1e.html
E_min, E_jmin = model_quartic_as_elliptic_curve(20, -9, 49, 318, 23, 106, -59, 81)
E_min.hyperelliptic_polynomials()[0].roots()
[(978559, 1)]
p0 = E_min(978559, 0)
(978559 : 0 : 1)
p0 in E_min
True
p2=[QQ(47971729)/49, QQ(16603172706)/343]
p2 in E_min
True
p2 = E_min(*p2) # Convert from list to point
p3=[1237921, 1244044242]
p3 in E_min
True
p3 = E_min(*p3) # Convert from list to point
E_min.gens(algorithm='pari', pari_effort=10)
[(37313463163849/15896169 : 246260871807590439226/63378025803 : 1), 
(1684027681/81 : 69276074469938/729 : 1), 
(101452934809/9 : 32314461714969994/27 : 1)]
g1, g2, g3 = E_min.gens()
p0-g2 == p3
True
p0-g3==p2
True
"""

def k_to_abcd(m, n, k_num, k_den, x0_num, x0_den, y0_num, y0_den):
    # 1. Get the parameterized x and y from the conic
    x_sym, y_sym = parameterize_conic_solution(m, n, x0_num, x0_den, y0_num, y0_den)
    
    # 2. Evaluate x and y at the specific k
    k_val = QQ(k_num) / k_den
    # We use subs(k=...) to get the rational values
    k = var('k')
    xv = x_sym.subs(k=k_val)
    yv = y_sym.subs(k=k_val)
    
    # 3. Calculate t from the second conic (Eq 4)
    a0, a2, a4, a6 = mn_to_coeffs_t2(m, n)
    t_sq = (a2*xv**2 + a4*xv + a6) / a0
    tv = t_sq.sqrt()
    
    # 4. Convert to the Elkies/Tomita variables r, s, t
    # r = x + y, s = x - y
    rv = xv + yv
    sv = xv - yv
    
    # 5. Extract A, B, C, D
    # These are the numerators after finding a common denominator
    all_fracs = [rv, sv, tv, QQ(1)]
    common_den = lcm([f.denominator() for f in all_fracs])
    
    A = abs(rv * common_den)
    B = abs(sv * common_den)
    C = abs(tv * common_den)
    D = abs(common_den)
    
    return sorted([Integer(A), Integer(B), Integer(C)]), Integer(D)
"""
% python -m elliptical.curves_tomita k_to_abcd 20 -9 -59 81 49 318 23 106
([95800, 217519, 414560], 422481)
"""

def get_t_from_k_generalized(k, x_val, m, n):
    M = m*m + m*n + n*n
    N = m*m - n*n
    L = 2*m*n + n*n
    
    # Plug the derived x into the second quadratic for t:
    # Mt^2 = Lx^2 - 2Nx - M
    t_sq = (L*x_val*x_val - 2*N*x_val - M) / M
    
    # We return the rational square root. 
    # Because we are solving A^4 + B^4 + C^4 = D^4, the sign of t 
    # doesn't change the final solution.
    return sqrt(t_sq)
 
def get_parametrized_xy(k, m, n):
    M = m*m + m*n + n*n
    N = m*m - n*n
    L = 2*m*n + n*n
    
    # We use the point (1, 0) to parameterize the first quadratic.
    # The line through (1, 0) with slope k is: y = k(x - 1)
    # Substituting this into My^2 = -Nx^2 - Lx + M and solving for x:
    
    x_val = (M*k*k - L - N) / (M*k*k + N)
    y_val = k * (x_val - 1)
    
    return x_val, y_val

def generate_euler_solution(E_min, P, m, n, k0, y0_quartic):
    a1, a2, a3, _, _ = E_min.ainvs()
    
    # Check for Point at Infinity
    if P.is_zero():
        k = k0 
    else:
        denom = (2 * P[1] - a1 * P[0] - a3)
        if denom == 0:
            # This is the point that maps to the 'infinity' of the parameter u
            # For these quartics, that usually means k is effectively infinite.
            # You can handle this by taking the limit in your x, y formulas.
            return "k is infinite at this point (Torsion Point)"
        
        u = (2 * y0_quartic * (P[0] + a2)) / denom
        k = u + k0
    print('k', k)

    # 2. Derive Parametrization Coefficients from (m, n)
    # These are the values from Elkies' original derivation
    M = m*m + m*n + n*n
    N = m*m - n*n
    L = 2*m*n + n*n
    
    # Tomita's x, y parametrization (Equation 5 & 6 generalized)
    # These coefficients depend on (m, n)
    q_denom = (M + L)*k*k + M + N
    
    # These numerators are the results of the parametrization of the 
    # specific quadratic: (m^2+mn+n^2)y^2 = -(m^2-n^2)x^2 ... 
    # To be fully general, these should be pulled from your 
    # derive_tomita_quartic function.
    
    # For now, we calculate r, s, t directly from the k-mapping:
    # Based on Elkies (1988), r, s, t are rational functions of k
    # that satisfy r^4 + s^4 + t^4 = 1.
    
    tx, ty = get_parametrized_xy(k, m, n)
    
    r_rat = tx + ty
    s_rat = tx - ty
    
    # t is the third variable satisfying the fourth-power relation
    # r^4 + s^4 + t^4 = 1 => t = (1 - r^4 - s^4)^(1/4)
    # In practice, we solve the second quadratic from Tomita's method:
    t_rat = get_t_from_k_generalized(k, tx, m, n)

    # 3. Clear Denominators
    common_denom = lcm([r_rat.denominator(), s_rat.denominator(), t_rat.denominator()])
    
    A = abs(r_rat * common_denom)
    B = abs(s_rat * common_denom)
    C = abs(t_rat * common_denom)
    D = abs(common_denom)
    
    return sorted([int(A), int(B), int(C)]), int(D)
"""
a4, a6 = 2265722465761, -3154189403034549278
E_min = EllipticCurve([a4, a6])
E_min.hyperelliptic_polynomials()[0].roots()
[(978559, 1)]
p0 = E_min(E_min.hyperelliptic_polynomials()[0].roots()[0][0], 0)
E_min.gens(algorithm='pari', pari_effort=10)
[(37313463163849/15896169 : 246260871807590439226/63378025803 : 1), (1684027681/81 : 69276074469938/729 : 1), (101452934809/9 : 32314461714969994/27 : 1)]
g1, g2, g3 = E_min.gens()
p3=p0-g1
p2=p0-g3
k0 = QQ(-59)/81
y0 = QQ(1222021742)/6561
generate_euler_solution(E_min, p0, 20, -9, k0, y0)
*** ZeroDivisionError: rational division by zero
>>> p0+g1
(98822808873914209/20881117009 : -32155995292015599252526194/3017384051151527 : 1)
generate_euler_solution(E_min, p0+g1, 20, -9, k0, y0)
k -1257455645286204321843837867740/15069677507922453335416025631
([469127667922560098297967958598520246471999154221606609719286700, 482733373780035471202275436545865297707320586482923557454835620, 957462603693495127512850795460242435562889242770658300285114836], 476012048020986994164337616826228853342372800193323421559450959)
p0+g2
(1237921 : -1244044242 : 1)
>>> generate_euler_solution(E_min, p0+g2, 20, -9, k0, y0)
k -759355832155450/4081087135881
([24635685568609618074908577397930, 24953554150672912060729156775230, 49877590627934529071542488767370], 24795474037173211434214437796837)
p0+g3
(47971729/49 : -16603172706/343 : 1)
>>> generate_euler_solution(E_min, p0+g3, 20, -9, k0, y0)
k -205218408718842700/54466708062033
([12672480311695159283691469649940499540, 12680505807785629799896054924632633340, 25499962273863820101935509553733202659], 12676494124757795615982047176103985391)

"""

def search_small_solutions_4(E_min, m, n, k0, y0):
    p0 = E_min(E_min.hyperelliptic_polynomials()[0].roots()[0][0], 0)
    gs = E_min.gens(algorithm='pari', pari_effort=10)
    g_set = set(gs)
    for g1 in gs:
        for g2 in gs:
            g_set.update([g1 + g2, g1 - g2])
    p_set = g_set | {p0 + g for g in g_set}
    p_list = sorted(p_set)
    print(f'Searching {len(p_list)} points')
    for p in p_list:
        abc, d = generate_euler_solution(E_min, p, m, n, k0, y0)
        if d < 1e27:
            print(f"Found Solution at point {p}:")
            print(abc, d)
"""
search_small_solutions_4(E_min, 20, -9, k0, y0)
Searching 32 points
k -59/81
Found Solution at point (0 : 1 : 0):
[930320, 2501002, 3957781] 3140740
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "<stdin>", line 12, in search_small_solutions_4
ValueError: too many values to unpack (expected 2)
>>> 930320**4 + 2501002**4 +  3957781**4 == 3140740**4
False
"""


def DEBUG():
    import pdb; pdb.set_trace()
    pass; pass; pass # opportumity to debug

"""
Follow section 2 of http://www.maroon.dti.ne.jp/fermat/dioph4e.html
for MacCloud's solution:
6306626244+ 2751562404+ 2190764654=6385232494

mn_to_coeffs_y2(8, -5)
(153, -779, -206, 80)
 mn_to_coeffs_t2(8, -5)
(153, 412, -320, -103)

COMTINUE $$$$$$$$$$$$$$$$$$$$$$$$$$
"""

if __name__ == "__main__":
    import sys
    command = eval(sys.argv[1])
    print(command)
    args = list(map(int, sys.argv[2:]))
    result = command(*args)
    print(result)

""" TODO
Why is search_small_solutions_4 broken?
"""