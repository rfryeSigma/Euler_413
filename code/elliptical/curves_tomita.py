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


#def model_quartic_as_elliptic_curve(m, n, x0_num, x0_den, y0_num, y0_den, pt_num, pt_den):
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
    
    E = EllipticCurve([a1, a2, a3, a4, a6])
    E_short = E.short_weierstrass_model()
    E_min = E_short.minimal_model()
    E_jmin = EllipticCurve(j=E.j_invariant()) # often shorter, but quadratic twist
    return E_min, E_jmin, shifted_poly
"""
q_res = make_quartic(QQ(20/-9), (QQ(49/318), QQ(23/106)))
model_quartic_as_elliptic_curve(q_res[2], QQ(-59/81))
(Elliptic Curve defined by y^2 = x^3 + 2265722465761*x - 3154189403034549278 over Rational Field, 
 Elliptic Curve defined by y^2 = x^3 + 2265722465761*x - 3154189403034549278 over Rational Field)

which matches C1:X^3+2265722465761X-3154189403034549278=Y^2
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

# broken
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

# broken
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

# broken
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

""" TODO
Clean up code to use solutions_curves.
Compare EC at end of curves_tomita.md with model_quartic_as_elliptic_curve
Why is search_small_solutions_4 broken?
"""

if __name__ == "__main__":
    import sys
    command = eval(sys.argv[1])
    print(command)
    args = list(map(int, sys.argv[2:]))
    result = command(*args)
    print(result)

