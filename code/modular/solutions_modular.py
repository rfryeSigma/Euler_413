"""
Explore whether modular curves can help find new solutions to a^4 + b^4 + c^4 = d^4.
"""

from solutions import known
from sage.all import help, oo, pari, sage_eval, solve, var, \
    QQ, Rational, \
    Curve, DiagonalQuadraticForm, Jacobian,PolynomialRing
from sage.rings.polynomial.polynomial_rational_flint import Polynomial_rational_flint

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
    m, n = u.numerator(), u.denominator()
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

So give up on Elkies rules. Just restrict to odd and even.
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
            #print(f"Error: {e}")
            return None
        assert point is not None
            
        X, V, Z = point
        # Convert back to original x, y
        x_val = (X / Z - b) / (2 * a)
        y_val = V / Z
        assert a0 * y_val**2 == a * x_val**2 + b * x_val + c
        # print(x_val, y_val)
    return (x_val, y_val)

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

def scan_known(n_known: int=6, max_h: int=100_000_000,
               max_d: int=int(1e27)) -> None:
    """For all abcd in first n_known known solutions, generate 12 u.
    For each u with height less than max_h:
        Verify no obstructions on u
    
    TODO SEE NOTES AT END OF FILE
    """
    for val in list(known.values())[:n_known]:
        abcd = val['abcd']
        print(abcd[-1])
        for h, u in abcd_to_h_u(abcd):
            if h > max_h: continue
            #print(u, check_elkies_rules(u))
            assert find_point_instantly(QQ(1)/u) is not None, f'No point for u {u}'
            print('\t', u)
            #TODO generate JD2 (smarter code) for each u
            # I think I should use u, but maybe 1/u?
            # Decide by running on C0, C1, C5, C7, C9, Cx
    pass
"""
Use pieza_in_pairs.py
    u_to_D_to_EC(u)

Use solutions_curves.py
    make_quartic(mn: Rational, quad_xy: tuple)

Use curves_tomita.py
    model_quartic_as_elliptic_curve(quartic_poly, quartic_x: Rational)
Use modular.log
# The Sage-native way to convert a quartic to an elliptic curve
R.<x,y> = QQ[]
f = x^4 + x + 1  # your quartic
C = Curve(y^2 - f)
E = Jacobian(C)
print(E.minimal_model())

def model_quartic_to_elliptic(quartic_poly):
    # pari.ellfromquartic handles the transformation formulas internally
    # It expects a polynomial in one variable
    res = pari.ellfromquartic(quartic_poly)
    # Convert PARI object back to a Sage Elliptic Curve
    return EllipticCurve(res.ell_to_sage())

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
TODO
Scan first 6 known solutions.
Identify 12 u.
Verify no obstructions on u
Generate JD2 (smarter code) for each u
Search for square quartic point.
Get conductor factored, torsion points, gens
height_pairing_matrix, new points.
Check whether known solutions.

map u to E, conductor, number of solutions.
Check whether isogenous with another curve.

>>> abcd_to_small_u((5_870_000, 8_282_543, 11_289_040, 12_197_457), 9e99)
[-400/37, -1010819791893/158785763120, -338000022077/104832225800, -2433/920, 
-35798568240/28609248113, -93/80, -84237/359800, 20632147/117135680, 
8685847/22963880, 450668400/124346123, 11502160/2925527, 1867333/457280]
>>> JD2 = u_to_D_to_EC(QQ(-93/80))
>>> JD2[0]
Elliptic Curve defined by y^2 = x^3 + 1687750917943790881*x - 294299265667029867546450078 over Rational Field
>>> JD2[1]
876967441/40960000*v^4 - 1160148721/10240000*v^3 + 5260917/32000*v^2 - 930349361/5120000*v - 540248559/10240000
>>> JD2[0].conductor()
31173291851033505900611518136
>>> from sage.all import factor
>>> factor(JD2[0].conductor())
2^3 * 17 * 41 * 73 * 89 * 241 * 1249 * 2137 * 2887 * 6569 * 70537
>>> JD2[0].torsion_points()
[(0 : 1 : 0), (171390639 : 0 : 1)]
>>> JD2[0].hyperelliptic_polynomials()[0].roots()[0][0]
171390639
>>> pari(JD2[0]).ellrank(15)
[2, 2, 2, [[18821840650969/8649, 94094386506723420026/804357], [2216934161202292017321/77272324441, 104489172057583918416101354132502/21480083475784739]]]
ellrank returns different values each time.
[2, 2, 2, [[1057203921, 51687652673778], [7086871502336635716361489/30329632086478209, 56104185538669647074758265030211769294/5282028171881234074661823]]]
 
>> g = JD2[0].gens(algorithm="pari", pari_effort=15)
[(7086871502336635716361489/30329632086478209 : 56104185538669647074758265030211769294/5282028171881234074661823 : 1), (18821840650969/8649 : 94094386506723420026/804357 : 1)]
but sage cache the values.

[r,R,s,V] = ellrank(E)
W = ellsaturation(E, V, 100)
takes too long when V is longer than the rank.

How to do LLL basis reduction on the generators to get smaller height points.
# 1. Get the height matrix from your rank-2 (or rank-k) curve
M = E.height_pairing_matrix(gens)

# 2. Call PARI's qflll on that matrix
# We convert the Sage matrix to a PARI object first
pari_M = pari(M)
T = pari_M.qflllgram() 

# 3. Apply the transformation T to your generators
# If gens = [P1, P2], and T = [[t11, t12], [t21, t22]]
# New P1 = t11*P1 + t21*P2, etc.
new_gens = []
for col in range(len(gens)):
    # Combine original generators using the coefficients in the T matrix
    new_P = sum(int(T[i, col]) * gens[i] for i in range(len(gens)))
    new_gens.append(new_P)

# 4. new_gens now contains the LLL-reduced basis (shortest heights)

"""