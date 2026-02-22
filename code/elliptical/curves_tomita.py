"""
Code to explore how to manipulate elliptic curves like Seiji Tomita does.
See my notes in curves_tomita.md
and Tomita's notes in http://www.maroon.dti.ne.jp/fermat/dioph4e.html
"""

from sage.all import Integer, QQ, PolynomialRing, Rational, Conic, \
        var, solve, lcm, numerator, denominator

# Gemini claims that this converts x^4+y^4+z^4=1 to elliptical curve (X, Y)
# But I doubt that it is True
def abcd_to_XY(abcd: tuple) -> set:
    """Convert (a, b, c, d) representing a^4 + b^4 + c^4 = d^4 
    to 12 pairs of rational points (X, Y)
   """
    XY_set = set()
    A, B, C, D = map(QQ, abcd)
    xx, yy, zz = A/D, B/D, C/D
    for x, y, z in ((xx, yy, zz), (zz, xx, yy), (zz, yy, xx),):
        for sx in (x, -x):
            for sy in (y, -y):
                for sz in (z, -z):
                    xx_yy = (sx*sx + sy*sy)
                    X = 2 * (1 + sz*sz) / xx_yy
                    Y = 4 * sz * (sx*sx - sy*sy) / (xx_yy * xx_yy)
                    XY_set.add((X, Y))
    return XY_set

"""
abcd_to_XY((95_800, 217_519, 414_560, 422_481))
{
(112902355361/45259408400, 9343828200105588969/5121035121794976400), 
(375335670722/219174508961, 20163255227132615488800/48037465378295469299521), 
(700700377922/56492155361, 26717737392167136531840/3191363617331361040321), 
(375335670722/219174508961, -20163255227132615488800/48037465378295469299521), 
(112902355361/45259408400, -9343828200105588969/5121035121794976400), 
(700700377922/56492155361, -26717737392167136531840/3191363617331361040321)
}
"""

# Tomita's eq 1, 2 are Elkies eq 3b,c for y^2 and t^2
def mn_to_coeffs_y2(m: int, n: int) -> tuple:
    """Apply m, n to the formula
    (2m^2+n^2)y^2=-(6m^2-8mn+3n^2)x^2-2(2m^2-n^2)x-2mn
    and return the evaluated coefficients (a0, a2, a4, a6)
    """
    m2 = m*m
    n2 = n*n
    mn = m*n
    a0 = 2 * m2 + n2
    a2 = -(6 * m2 - 8 * mn + 3 *n2)
    a4 = -2 * (2 * m2 - n2)
    a6 = -2 * mn
    return (a0, a2, a4, a6)

def mn_to_coeffs_t2(m: int, n: int) -> tuple:
    """Apply m, n to the formula
    (2m2+n2)t2=4(2m2-n2)x2+8mnx+(n2-2m2)
    and return the evaluated coefficients (a0, a2, a4, a6)
    """
    m2 = m*m
    n2 = n*n
    mn = m*n
    a0 = 2 * m2 + n2
    a2 = 4 * (2 * m2 - n2)
    a4 = 8 * mn
    a6 = n2 - 2 * m2
    return (a0, a2, a4, a6)

""""
% python -m elliptical.curves_tomita mn_to_coeffs_y2 20 -9
(881, -4083, -1438, 360)

% python -m elliptical.curves_tomita mn_to_coeffs_t2 20 -9
(881, 2876, -1440, -719)
"""

def find_rational_points_in_conic(m: int, n: int) -> list:
    """Find rational points in Tomita eq 3.
    """
    R = PolynomialRing(QQ, 'x, y, z')
    x, y, z = R.gens()
    a0, a2, a4, a6 = mn_to_coeffs_y2(m, n)
    print(a0, a2, a4, a6)
    f = a0*y**2 -a2*x**2 -a4*x*z -a6*z**2
    C = Conic(f)
    pts = C.rational_points(bound=400)
    return pts

"""
% python -m elliptical.curves_tomita find_rational_points_in_conic 20 -9
<function find_rational_points_in_conic at 0x158fae8e0>
881 -4083 -1438 360
[(1/354 : -75/118 : 1), (1/354 : 75/118 : 1), (49/318 : -23/106 : 1), (49/318 : 23/106 : 1)]
"""

def parameterize_conic_solution(m: int, n: int, x0_num :int, x0_den: int,
                                y0_num: int, y0_den: int):
    """Paramterize the conic solution
    """
    R = PolynomialRing(QQ, 'x, y, z')
    x, y, z = R.gens()
    a0, a2, a4, a6 = mn_to_coeffs_y2(m, n)
    x0 = QQ(x0_num) / QQ(x0_den)
    y0 = QQ(y0_num) / QQ(y0_den)
    var('x_var, k')
    line = k*(x_var - x0) + y0
    eq = a0*(line)**2 -a2*x_var**2 + a4*x_var -a6
    roots = solve(eq == 0, x_var)
    x_k = roots[1].rhs().simplify_full()
    return x_k

"""
 % python -m elliptical.curves_tomita parameterize_conic_solution 20 -9 49 318 -23 106
<function parameterize_conic_solution at 0x161c0ea20>
1/318*(43169*k^2 + 60789*k + sqrt(43176288513*k^2 + 3474091350*k + 183791406681) + 228642)/(881*k^2 + 4083)

"""
def eval_xy_at_point(m: int, n: int, x_num: int, x_den: int, 
                     y_num: int, y_den: int,) -> bool:
    """Evaluate Tomita's eq 3 with m/n for (x, y) at given rational point
    """
    a0, a2, a4, a6 = mn_to_coeffs_y2(m, n)
    y = QQ(y_num) / QQ(y_den)
    lhs = a0 * y * y
    x = QQ(x_num) / QQ(x_den)
    rhs = (a2 * x + a4) * x + a6
    print(f'lhs: {lhs}')
    print(f'rhs: {rhs}')
    return lhs == rhs

"""
eval_xy_at_point(20, -9, 49, 318, 23, 106)
lhs: 466049/11236
rhs: 466049/11236
True
"""

if __name__ == "__main__":
    import sys
    command = eval(sys.argv[1])
    print(command)
    args = list(map(int, sys.argv[2:]))
    result = command(*args)
    print(result)

