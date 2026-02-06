"""
I am trying to understand mapping between abcd solutions 
and rational characteristics of their Elliptic Curves.

I am following Tito Piezas III's suggestions in
https://math.stackexchange.com/questions/4857229/on-why-solutions-to-x4y4z4-1-come-in-pairs
and in
https://math.stackexchange.com/questions/1853223/distribution-of-primitive-pythagorean-triples-ppt-and-of-solutions-of-a4b4

"""
from solutions import known
from itertools import combinations
from math import isqrt

def abcd_to_xyz(abcd: tuple) -> tuple:
    """Convert (a, b, c, d) representing a^4 + b^4 + c^4 = d^4 
    to rational form (x, y, z) representing x^4 + y^4 + z^4 = 1
    """
    A, B, C, D = map(QQ, abcd)
    x, y, z = A/D, B/D, C/D
    return x, y, z

def verify_uvw(u, v, w) -> bool:
    """Verify that (u, v, w) satisfy 2(u + v + w) - uvw - 4 = 0
    This is True for calculated triples and their permuations,
    but False for remixed triples.
    """
    return 2 * (u + v + w) - u * v * w - 4 == 0

def xyz_to_uvw(xyz: tuple) -> tuple:
    """Convert (x, y, z) representing x^4 + y^4 + z^4 = 1
    to EC parameters (u, v, w)
    """
    x, y, z = xyz
    x_y, y_z, z_x = x - y, y - z, z - x
    u = (x_y*x_y - z*z -1) / (x*x - x*y + y*y + x_y)
    v = (y_z*y_z - x*x -1) / (y*y - y*z + z*z + y_z)
    w = (z_x*z_x - y*y -1) / (z*z - z*x + x*x + z_x)
    assert verify_uvw(u, v, w), f'{u}, {v}, {w}'
    return u, v, w

# Seems to work best when a, b, c in increasing size
xyz_to_uvw(abcd_to_xyz((95_800, 217_519, 414_560, 422_481)))
(1000/47, -1041/320, -9/20)
xyz_to_uvw(abcd_to_xyz((414560, 217519, 95800, 422_481)))
(-71490240/101943281, -167767337/43538900, -6899820729/369596780)
xyz_to_uvw(abcd_to_xyz((-414560, -217519, -95800, 422_481)))
(-1041/320, 1000/47, -9/20)

def abcd_to_small_u(abcd: tuple, thresh: int=10_000) -> list:
    """Convert (a, b, c, d) to (x, y, z) and then consider
    6 permutations and 8 combinations of signs on (x, y, z).
    Since the u, v, w formulas rotate (x, y, z), only swap y, z.
    to generate the set of all EC parameters u.
    Select those with abs of denominator or numerator < threshold.
    """
    u_set = set()
    x, yy, zz = abcd_to_xyz(abcd)
    for y, z in ((yy, zz), (zz, yy)):
        for sx in (x, -x):
            for sy in (y, -y):
                for sz in (z, -z):
                    uvw = xyz_to_uvw((sx, sy, sz))
                    u_set.update(uvw)
    u_list = []
    for u in u_set:
        if abs(numerator(u)) < thresh or abs(denominator(u)) < thresh:
            u_list.append(u)
    return sorted(u_list)

"""The 11 solutions on curve C1:
1. 422481; 414560, 217519, 95800 (Frye, 1988); C1, u = -9/20.
abcd_to_small_u((414560, 217519, 95800, 422_481))
[-1041/320, -4209/3500, -9/20, 30080/6007, 1000/47]

2. 2813001; 2767624, 1390400, 673865 (MacLeod 1997); C1, u = -9/20.
abcd_to_small_u((2767624, 1390400, 673865, 2813001))
[-1425/412, -9/20, 34225/6692, 5728/215]

17. 1679142729; 1670617271, 632671960, 50237800 (Tomita, 2006); C1, u = -9/20.
abcd_to_small_u((1670617271, 632671960, 50237800, 1679142729))
[-1041/320, -9/20, 1000/47]

25. 15434547801; 15355831360, 5821981400, 140976551 (Tomita, 2007); C1
abcd_to_small_u((15355831360, 5821981400, 140976551, 15434547801))
[-1425/412, -9/20, 5728/215]

42. 29999857938609; 27239791692640, 22495595284040, 7592431981391 (Tomita, 2006); C1
abcd_to_small_u((27239791692640, 22495595284040, 7592431981391, 29999857938609))
[-4209/3500, -9/20, 30080/6007]

58 101783028910511968041; 99569174129827461335, 21710111037730547416, 54488888702794271560 (Piezas, 2024); C1
abcd_to_small_u((99569174129827461335, 21710111037730547416, 54488888702794271560, 101783028910511968041))
[-9/20]
75. 120175486227071990769561; 30248376090268690676600, 118508989446504950664160, 56915898438422390129561 (Tomita, 2024); C1
abcd_to_small_u((30248376090268690676600, 118508989446504950664160, 56915898438422390129561, 120175486227071990769561))
[-9/20]
80. 1171867103503245199920081; 1165970778032514255823760, 440517744543240750721000, 59421842165791512201169 (Tomita, 2024); C1
abcd_to_small_u((1165970778032514255823760, 440517744543240750721000, 59421842165791512201169, 1171867103503245199920081))
[-9/20]
84. 6714012701109174954871521; 1758067984180618846616200, 6632467268281371571709360, 3057432874236989781768479 (Tomita, 2024); C1
abcd_to_small_u((1758067984180618846616200, 6632467268281371571709360, 3057432874236989781768479, 6714012701109174954871521))
[-9/20]
87. 21291952935426564624339201; 5328636655728999148343576, 20991236668646283695879935, 10137374115207940432133560 (Tomita, 2024); C1
abcd_to_small_u((5328636655728999148343576, 20991236668646283695879935, 10137374115207940432133560, 21291952935426564624339201))
[-9/20]
92. 227529118288906398066378489, 85818832944459457142858489, 226369052354324181334408840, 2650718685573298353948640 (Piezas, 2024); C1
abcd_to_small_u((85818832944459457142858489, 226369052354324181334408840, 2650718685573298353948640, 227529118288906398066378489))
[-9/20]
"""

def is_square(u) -> bool:
    """Return whether a rational is a perfect square"""
    n = numerator(u)
    if isqrt(n)**2 != n: return False
    d = denominator(u)
    return isqrt(d)**2 == d

def un_square(u):
    """Return isqrt of numerator and denominator of a rational"""
    n = isqrt(numerator(u))
    d = isqrt(denominator(u))
    r = QQ(n) / QQ(d)
    assert r * r == u, f'{u}, {r}'
    return r

def uv_to_Pk(u, v):
    """Calculate the EC parameters Pk from u, v
    """
    u2 = u * u
    u3 = u2 * u
    v2 = v * v
    v3 = v2 * v
    P0 =(2 + u2) * (2 + v2) * (12 - 8*u + 2*u2 - 8*v + 2*v2 + u2*v2)
    P1 = (-4 + 4*u + 2*v - u2*v) * (8*u -4*u2 + 4*v -8*u*v + 2*u2*v -4*v2 + 2*v3 + u2*v3)
    P2 = 2 * (-4 + 4*u + 2*v -u2*v) * (4 -2*u -4*v +u*v2) * (-u + v)
    P3 = (4 - 2*u - 4*v + u*v2) * (4*u - 4*u2 + 2*u3 + 8*v - 8*u*v - 4*v2 + 2*u*v2 + u3*v2)
    return (P0, P1, P2, P3)

def uv_to_D(u, v):
    """Calculate a new EC parameter D from u, v
    """
    u2 = u * u
    u3 = u2 * u
    u4 = u2 * u2
    d0 = 4 * (-6 - 2 * u + u2) * (2 - 2 * u + u2) 
    d1 = - 8 * (-2 - 4 * u + u2) * (2 - 2 * u + u2)
    d2 = - 16 * u * (4 - 3 * u + u2)
    d3 = - 4 * (4 - 12 * u + 4 * u2 - 2 * u3 + u4) 
    d4 = (4 - 8 * u - 4 * u3 + u4)
    v2 = v * v
    v3 = v2 * v
    v4 = v2 * v2
    D2 = d0 + d1 * v + d2 * v2 + d3 * v3 + d4 * v4
    if is_square(D2):
        return un_square(D2)
    else:
        print(f'u {u}, v {v} to D2 {D2} is not a perfect square ')
    return None

u = QQ(-9) / QQ(20)
v = QQ(-1041) / QQ(320)
uv_to_Pk(u, v)
(3037974602413721121/1677721600000000,
 2903996886978203679/1677721600000000,
 36016220867741673/41943040000000,
 -2059498434165383679/1677721600000000)
uv_to_D(u, v)
2126704839/40960000

def uv_to_xyz(u, v):
    """Calculate (x, y, z) from u, v
    """
    P0, P1, P2, P3 = uv_to_Pk(u, v)
    D = uv_to_D(u, v)
    if D is None: return None
    u2 = u * u
    v2 = v * v
    denom = P0 + D*D
    x1 = (u2 - 2*u + 2) * (v2 - 2*v)
    y1 = 2*(u + v - 2) * (u + v - u*v)
    z1 = (v2 - 2*v + 2) * (u2 - 2*u)
    pair = [0]*2
    for inx, sD in enumerate((D, -D)):
        x = (P1 + x1 * sD) / denom
        y = (-(P1 + P2 + P3) + y1 * sD) / denom
        z = (P3 + z1 * sD) / denom
        pair[inx] = (x, y, z)
    return pair

u = QQ(-9) / QQ(20)
v = QQ(-1041) / QQ(320)
uv_to_xyz(u,v)
[(1670617271/1679142729, 632671960/1679142729, -50237800/1679142729),
 (-95800/422481, -414560/422481, -217519/422481)]

def abcd_to_small_xyz(abcd: tuple, thresh: int=10_000) -> list:
    """Convert (a, b, c, d) to small u.
    For all permutations of the small u taken 2 at a time,
    solve for the corresponding (x, y, z) pairs.
    """
    xyz_list = []
    for u, v in combinations(abcd_to_small_u(abcd), 2):
        xyz = uv_to_xyz(u, v)
        if xyz is None: continue
        for x, y, z in xyz:
            xyz_list.append((x, y, z))
    return xyz_list

abcd_to_small_xyz((414560, 217519, 95800, 422_481))
# u -1041/320, v -4209/3500 to D2 14290666362400669957976538801/1573519360000000000000000 is not a perfect square 
# u -1041/320, v 30080/6007 to D2 198909725101194840504881377761/3413268476026948157440000 is not a perfect square 
# u -4209/3500, v 1000/47 to D2 93428880408453190126721463/26152040359375000000 is not a perfect square 
# u 30080/6007, v 1000/47 to D2 11211379600882246264097552592/6353630573412954106081 is not a perfect square 
[(217519/422481, 414560/422481, 95800/422481),
 (50237800/1679142729, -632671960/1679142729, -1670617271/1679142729),
 (632671960/1679142729, -50237800/1679142729, 1670617271/1679142729),
 (-414560/422481, -217519/422481, -95800/422481),
 (27239791692640/29999857938609,
  22495595284040/29999857938609,
  -7592431981391/29999857938609),
 (217519/422481, -95800/422481, -414560/422481),
 (95800/422481, -217519/422481, 414560/422481),
 (-22495595284040/29999857938609,
  -27239791692640/29999857938609,
  7592431981391/29999857938609),
 (22495595284040/29999857938609,
  -7592431981391/29999857938609,
  27239791692640/29999857938609),
 (-95800/422481, -414560/422481, 217519/422481),
 (414560/422481, 95800/422481, 217519/422481),
 (-632671960/1679142729, -1670617271/1679142729, 50237800/1679142729)]
# The unique denominators in the result are
# 422481 (1 in table)
# 1679142729 (7 in table)
# 29999857938609 (12 in table)
