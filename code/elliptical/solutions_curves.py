"""
Use the uv and their D in solutions_uv.csv with Tomita's elliptic curves
to generates solution.
See my notes in curves_tomita.{md/py}
and Tomita's notes in http://www.maroon.dti.ne.jp/fermat/dioph4e.html
"""
#from sage.all import Integer, QQ, lcm, ceil, gcd, solve, sqrt, var, \
#        EllipticCurve, parallel, PolynomialRing
import csv

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
{(-71490240, 101943281), (-167767337, 43538900), (-6899820729, 369596780)}
"""

def DEBUG():
    import pdb; pdb.set_trace()
    pass; pass; pass # opportumity to debug

"""
"""

if __name__ == "__main__":
    import sys
    command = eval(sys.argv[1])
    print(command)
    args = list(map(int, sys.argv[2:]))
    result = command(*args)
    print(result)

""" TODO
Use u to generate solutions and correlate with map.

"""