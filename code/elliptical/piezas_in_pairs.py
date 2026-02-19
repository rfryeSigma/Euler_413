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
from sage.all import Integer, QQ, Rational, gcd, lcm, numerator, denominator

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
    to parameters (u, v, w)
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

"""The 12 solutions on curve C1:
1. 422481; 414560, 217519, 95800 (Frye, 1988); C1, u = -9/20.
[-1041/320, -4209/3500, -9/20, 30080/6007, 1000/47]

2. 2813001; 2767624, 1390400, 673865 (MacLeod 1997); C1, u = -9/20.
[-1425/412, -9/20, 34225/6692, 5728/215]

17. 1679142729; 1670617271, 632671960, 50237800 (Tomita, 2006); C1, u = -9/20.
[-1041/320, -9/20, 1000/47]

23. 15434547801; 15355831360, 5821981400, 140976551 (Tomita, 2007); C1
[-1425/412, -9/20, 5728/215]

42. 29999857938609; 27239791692640, 22495595284040, 7592431981391 (Tomita, 2006); C1
[-4209/3500, -9/20, 30080/6007]

46. 573646321871961; 514818101299289, 440804942580160, 130064300991400 (Tomita, 2008); C1
[-9/20, 34225/6692]

58 101783028910511968041; 99569174129827461335, 21710111037730547416, 54488888702794271560 (Piezas, 2024); C1
[-9/20]

75. 120175486227071990769561; 30248376090268690676600, 118508989446504950664160, 56915898438422390129561 (Tomita, 2024); C1
[-9/20]

80. 1171867103503245199920081; 1165970778032514255823760, 440517744543240750721000, 59421842165791512201169 (Tomita, 2024); C1
[-9/20]

84. 6714012701109174954871521; 1758067984180618846616200, 6632467268281371571709360, 3057432874236989781768479 (Tomita, 2024); C1
-9/20]

87. 21291952935426564624339201; 5328636655728999148343576, 20991236668646283695879935, 10137374115207940432133560 (Tomita, 2024); C1
[-9/20]

92. 227529118288906398066378489, 85818832944459457142858489, 226369052354324181334408840, 2650718685573298353948640 (Piezas, 2024); C1
[-9/20]
"""

def is_square(u) -> bool:
    """Return whether a rational is a perfect square"""
    n = numerator(u)
    if 0 >= n or isqrt(n)**2 != n: return False
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

def uv_to_D(u, v, verbose=False):
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
        if verbose:
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
        pair = uv_to_xyz(u, v)
        if pair is None: continue
        for x, y, z in pair:
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

def known_to_unknown(max_d: int=int(1e27)) -> None:
    """For all abcd in known, make new xyz from combinations of u.
    Check whether the denominators are below max_d and in known.
    Report any new xyz.
    """
    new_xyz = list()
    new_denoms = set()
    known_denoms = set()
    for val in known.values():
        abcd = val['abcd']
        known_denoms.add(abcd[-1])
    for val in known.values():
        abcd = val['abcd']
        for u, v in combinations(abcd_to_small_u(abcd, thresh=9e99), 2):
            pair = uv_to_xyz(u, v)
            if pair is None: continue
            for xyz in pair:
                d = lcm([denominator(x) for x in xyz])
                if d >= max_d or d in known_denoms: continue
                if d in new_denoms: continue
                new_denoms.add(d)
                new_xyz.append(xyz)
                print(f'Found new {xyz} from {abcd} with u {u}, v {v}')
    print(f'Found {len(new_xyz)} new xyz')
    return sorted(new_xyz)

# known_to_unknown()
# Found new (535607712470322407570378200/918310142223097801006397529, 250530677254015598463660440/918310142223097801006397529, -889106559932545369165762471/918310142223097801006397529) from (95800, 217519, 414560, 422481) with u -1041/320, v 52463/660460
# Found new (118194421251475239056505903/123188180833923372056627153, 66491673395168374249746120/123188180833923372056627153, -62831773759131557571594880/123188180833923372056627153) from (2164632, 31669120, 41084175, 44310257) with u -12065/12396, v -84558637/193874100
# Found new (46402888024739111034420161/174088703841632292189275073, -92419682545114696981174360/174088703841632292189275073, -170289556324371670328363560/174088703841632292189275073) from (27450160, 108644015, 146627384, 156646737) with u -136/133, v -63528125/85096232
# Found new (411177854471028470696556192/504068891841730072306483681, 106958136069417067994530335/504068891841730072306483681, -435117527990435060232042280/504068891841730072306483681) from (39110088360, 49796687200, 71826977313, 76973733409) with u -267904/221337, v -1245/5012
# Found 4 new xyz
[(46402888024739111034420161/174088703841632292189275073,
  -92419682545114696981174360/174088703841632292189275073,
  -170289556324371670328363560/174088703841632292189275073),
 (535607712470322407570378200/918310142223097801006397529,
  250530677254015598463660440/918310142223097801006397529,
  -889106559932545369165762471/918310142223097801006397529),
 (411177854471028470696556192/504068891841730072306483681,
  106958136069417067994530335/504068891841730072306483681,
  -435117527990435060232042280/504068891841730072306483681),
 (118194421251475239056505903/123188180833923372056627153,
  66491673395168374249746120/123188180833923372056627153,
  -62831773759131557571594880/123188180833923372056627153)]
"""
123_188_180_833_923_372_056_627_153;
118194421251475239056505903,
66491673395168374249746120,
62831773759131557571594880

174_088_703_841_632_292_189_275_073;
170289556324371670328363560,
92419682545114696981174360,
46402888024739111034420161

504_068_891_841_730_072_306_483_681;
435117527990435060232042280,
411177854471028470696556192, 
106958136069417067994530335

918_310_142_223_097_801_006_397_529;
889106559932545369165762471,
535607712470322407570378200,
250530677254015598463660440
"""

def given_to_unknown(given: tuple, max_d: int=int(1e27)) -> None:
    """For given abcd, make new xyz from combinations of u.
    Check whether the denominators are below max_d and in known.
    Report any new xyz.
    """
    a, b, c, d = given
    assert 0 < a < b < c < d
    assert gcd(given) == 1
    assert a**4 + b**4 + c**4 == d**4
    new_denoms = set((d,))
    new_xyz = list()
    known_denoms = set()
    for val in known.values():
        abcd = val['abcd']
        known_denoms.add(abcd[-1])
    for u, v in combinations(abcd_to_small_u(given, thresh=9e99), 2):
        pair = uv_to_xyz(u, v)
        if pair is None: continue
        for xyz in pair:
            d = lcm([denominator(x) for x in xyz])
            print(f'denom {d}')
            if d in known_denoms: 
                print('\tknown')
                continue
            if d >= max_d:
                print(f'\tbig {float(d):.4e}')
                continue
            if d in new_denoms: 
                print('\trepeat')
                continue
            new_denoms.add(d)
            new_xyz.append(xyz)
            print(f'Found new {xyz} from {abcd} with u {u}, v {v}')
    print(f'Found {len(new_xyz)} new xyz')
    return sorted(new_xyz)

def known_to_unknown_inv(max_d: int=int(1e27)) -> None:
    """For all abcd in known, make new xyz from combinations of 
    inverses of u.
    Check whether the denominators are below max_d and in known.
    Report any new xyz.
    """
    one = QQ(1)
    new_xyz = list()
    new_denoms = set()
    known_denoms = set()
    found_big = 0
    found_raw = 0
    for val in known.values():
        abcd = val['abcd']
        known_denoms.add(abcd[-1])
    for val in known.values():
        abcd = val['abcd']
        for uu, vv in combinations(abcd_to_small_u(abcd, thresh=9e99), 2):
            for u, v in ((uu, one/vv), (one/uu, vv), (one/uu, one/vv)):
                pair = uv_to_xyz(QQ(1)/u, QQ(1)/v)
                if pair is None: continue
                for xyz in pair:
                    d = lcm([denominator(x) for x in xyz])
                    if d >= max_d:
                        found_big += 1
                        continue
                    found_raw += 1
                    if d in known_denoms or d in new_denoms: continue
                    new_denoms.add(d)
                    new_xyz.append(xyz)
                    print(f'Found new {xyz} from {abcd} with u {u}, v {v}')
    print(f'Found {found_big} too big and {found_raw} ok'
          f', but only {len(new_xyz)} new xyz')
    return sorted(new_xyz)

# known_to_unknown_inv()
# Found 1962 too big and 2694 ok, but only 0 new xyz

def known_to_unknown_all_u(max_d: int=int(1e27)) -> None:
    """For all abcd in known, collect u. From combinations of all u.
    Check whether the denominators are below max_d and in known.
    Report any new xyz.
    """
    new_xyz = list()
    new_denoms = set()
    known_denoms = set()
    u_set = set()
    found_raw = 0
    found_big = 0
    for val in known.values():
        abcd = val['abcd']
        known_denoms.add(abcd[-1])
        for u in abcd_to_small_u(abcd, thresh=9e99):
            u_set.add(u)
    for u, v in combinations(u_set, 2):
        pair = uv_to_xyz(u, v)
        if pair is None: continue
        for xyz in pair:
            d = lcm([denominator(x) for x in xyz])
            if d >= max_d: 
                found_big += 1
                continue
            found_raw += 1
            if d in known_denoms or d in new_denoms: continue
            new_denoms.add(d)
            new_xyz.append(xyz)
            print(f'Found new {xyz} with u {u}, v {v}')
    print(f'Found {found_big} too big and {found_raw} ok'
          f', but only {len(new_xyz)} new xyz')
    return sorted(new_xyz)

#known_to_unknown_all_u()
#Found 1962 too big and 2328 ok, but only 0 new xyz

def check_new_d_size() -> None:
    """For all abcd in known, make new xyz from combinations of u.
    Check whether any new denominators are below current d.
    """
    count = 0
    for i, val in enumerate(known.values(), start=1):
        abcd = val['abcd']
        d = abcd[-1]
        for u, v in combinations(abcd_to_small_u(abcd, thresh=9e99), 2):
            pair = uv_to_xyz(u, v)
            if pair is None: continue
            for xyz in pair:
                new_d = lcm([denominator(x) for x in xyz])
                if new_d >= d: continue
                count += 1
    return count

#check_new_d_size()
#183
# This means that solutions larger than known can find solutions less than 1e27.
# But it's not clear that any will be new, since they were derived from known.

def large_known_to_unknown(min_d: int=int(1e27), max_d: int=int(1e40)) -> None:
    """For all abcd in known, make new xyz from combinations of u.
    Keep those in range min_d <= new_d < max_d
    Use these to find new xvz and report those with new_d < min_d,
    but not in known.
    """
    known_denoms = set()
    big_denoms = set()
    big_xyz = list()
    new_denoms = set()
    new_xyz = list()
    # Collect big xyz from known abcd, and their denominators.
    for val in known.values():
        abcd = val['abcd']
        known_denoms.add(abcd[-1])
        for u, v in combinations(abcd_to_small_u(abcd, thresh=9e99), 2):
            pair = uv_to_xyz(u, v)
            if pair is None: continue
            for xyz in pair:
                d = lcm([denominator(x) for x in xyz])
                if d < min_d or d >= max_d: continue
                if d in big_denoms: continue
                big_denoms.add(d)
                big_xyz.append(xyz)
    print(f'Found {len(big_denoms)} big xyz') 
    # Collect more big xyz from big_xyz  
    for xyz in big_xyz:
        d = lcm([denominator(x) for x in xyz])
        assert all(d == denominator(x) for x in xyz)
        abcd = [numerator(x) for x in xyz] + [d]
        for u, v in combinations(abcd_to_small_u(abcd, thresh=9e99), 2):
            pair = uv_to_xyz(u, v)
            if pair is None: continue
            for xyz in pair:
                d = lcm([denominator(x) for x in xyz])
                if d < min_d or d >= max_d: continue
                if d in big_denoms: continue
                big_denoms.add(d)
                big_xyz.append(xyz)
    print(f'Found {len(big_denoms)} big xyz') 
    # Collect new small xyz from all big_xyz  
    for xyz in big_xyz:
        d = lcm([denominator(x) for x in xyz])
        assert all(d == denominator(x) for x in xyz)
        abcd = [numerator(x) for x in xyz] + [d]
        for u, v in combinations(abcd_to_small_u(abcd, thresh=9e99), 2):
            pair = uv_to_xyz(u, v)
            if pair is None: continue
            for xyz in pair:
                d = lcm([denominator(x) for x in xyz])
                if d >= min_d or d in known_denoms: continue
                if d in new_denoms: continue
                new_denoms.add(d)
                new_xyz.append(xyz)
                print(f'Found new {xyz} from {abcd} with u {u}, v {v}')
    print(f'Found {len(new_xyz)} new xyz')
    return sorted(new_xyz)

#large_known_to_unknown()
#Found 77 big xyz
#Found 78 big xyz
#Found 0 new xyz

def known_to_unknown_all_uu(max_d: int=int(1e27)) -> None:
    """For all abcd in known, collect u. For all (u, u),
    check whether the denominators are below max_d and in known.
    Report any new xyz.
    """
    new_xyz = list()
    new_denoms = set()
    known_denoms = set()
    u_set = set()
    found_raw = 0
    found_big = 0
    found_known = 0
    for val in known.values():
        abcd = val['abcd']
        known_denoms.add(abcd[-1])
        abcd = val['abcd']
        for u in abcd_to_small_u(abcd, thresh=9e99):
            u_set.add(u)
    for u in u_set:
        pair = uv_to_xyz(u, u)
        if pair is None: continue
        for xyz in pair:
            d = lcm([denominator(x) for x in xyz])
            if d >= max_d: 
                found_big += 1
                continue
            found_raw += 1
            if d in new_denoms: continue
            if d in known_denoms: 
                found_known += 1
                continue
            new_denoms.add(d)
            new_xyz.append(xyz)
            print(f'Found new {xyz} with u {u}, v {u}')
    print(f'Found {found_big} too big and {found_raw} ok'
          f' and {found_known} in known'
          f', but only {len(new_xyz)} new xyz')
    return sorted(new_xyz)

#known_to_unknown_all_uu()
#Found 0 too big and 0 ok and 0 in known, but only 0 new xyz
