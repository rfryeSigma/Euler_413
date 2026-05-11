""" Search for a recurrence pattern in (u, v).
1. Scatter plot from a suggestion by Phil Hadley

2. Lattice Dependency Check
    Final note: Why I gave up

    To check for a hidden "lattice structure" among 100 solutions, 
    we can use the **LLL algorithm** to see if the numerators and 
    denominators of $u$ and $v$ values satisfy a small linear relation. 
    If they do, they are traveling along specific geometric "tracks" 
    (rational curves) on the K3 surface.
"""
import csv
from pdb import runcall, set_trace
from sage.all import QQ, ZZ, Matrix, sage_eval, scatter_plot

""" ----------------
2. Scatter Plot
"""
def scatter_known(row_s :int=1, row_e :int=100, clip :float=50.0,
        file_name: str='solutions_uv.csv') -> None:
    """ Check known (u, v) pairs for recurrence patterns
    For initial data row after header use row_s = 2 .
    Each (u, v) pair covers 2 rows, so skips odd row numbers.
    """
    data = []
    with open(file_name, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file)
        for _ in range(1, row_s):
            next(reader)
        for i_row in range(row_s, row_e + 1):
            try:
                row = next(reader)
            except StopIteration: break
            if i_row%2 == 1: continue # duplicate (u, v)
            u_n, u_d, v_n, v_d = [int(row[i]) for i in range(4)]
            if not (-clip < float(u_n)/u_d < clip): continue
            if not (-clip < float(v_n)/v_d < clip): continue
            data.append((QQ(u_n)/u_d, QQ(v_n)/v_d))
    p = scatter_plot(data, marker='o', markersize=2, facecolor='black',
                 title='(u, v) scatter')
    p.show(gridlines=True)
    p.save('uv_scatter.png')

""" ----------------
2. Lattice Dependency Check
"""
def check_uv_lattice(p: int, q: int, r :int, s : int) -> tuple:
    """
    Checks if the components of (u, v) satisfy a small linear relation.
    u = p/q, v = r/s
    """
    # We construct a lattice with the components.
    # The large constant K forces the algorithm to find a relation that
    # sums exactly (or nearly) to zero.
    K = 10**40

    """
    # Trivial Linear: ap + bq + cr + ds = 0
    # Allows qp - pq = 0
    M = Matrix(ZZ, [
        [1, 0, 0, 0, p * K],
        [0, 1, 0, 0, q * K],
        [0, 0, 1, 0, r * K],
        [0, 0, 0, 1, s * K],
    ])
    """

    """
    # Linear Fractional Mobius Transform: v = (au+b) / (cu+d)
    # Allows -qp + pq = 0
    M = Matrix(ZZ, [
        [1, 0, 0, 0, p * s * K],
        [0, 1, 0, 0, r * q * K],
        [0, 0, 1, 0, q * s * K],
        [0, 0, 0, 1, p * r * K],
    ])
    """

    # Quadratic Relation: 
    # Can still get trivial solutions when height of v much larger than u.
    M = Matrix(ZZ, [
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, p * 1 * K],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, p * p * K],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q * 1 * K],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q * q * K],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, p * q * K],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, r * 1 * K],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, r * r * K],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, s * 1 * K],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, s * s * K],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, r * s * K],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, p * r * K],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, p * s * K],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, q * r * K],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, q * s * K],
    ])
    
    # The first row of the LLL-reduced matrix is the shortest vector
    L = M.LLL()
    for i in range(L.nrows()):
        relation = L[0][:-1] # remove the 'K' column
        residual = L[0][-1]  # should be 0 if a perfect relation exists
        if 0 != residual: break # no good relation
        if not (2 < relation.hamming_weight() < 6):
            residual = -1
            break # too simple or too complex
        u_any = 0 < relation[:5].hamming_weight()
        v_any = 0 < relation[5:10].hamming_weight()
        inter = 0 < relation[10:].hamming_weight()
        if not (inter or (u_any and v_any)):
            residual = -1
            break # does not include both u and v
    return relation, residual

def check_known(row_s :int=1, row_e :int=100,
        file_name: str='solutions_uv.csv') -> None:
    """ Check known (u, v) pairs for recurrence patterns
    For initial data row after header use row_s = 2 .
    Each (u, v) pair covers 2 rows, so skips odd row numbers.
    """
    with open(file_name, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file)
        for _ in range(1, row_s):
            next(reader)
        for i_row in range(row_s, row_e + 1):
            try:
                row = next(reader)
            except StopIteration: break
            if i_row%2 == 1: continue # duplicate (u, v)
            u_n, u_d, v_n, v_d = [int(row[i]) for i in range(4)]
            rel, res = check_uv_lattice(u_n, u_d, v_n, v_d)
            if res != 0: continue
            print(u_n, u_d, v_n, v_d, '->', rel)
"""
python -um modular.uv_recurrence check_known 20
<function check_known at 0x102a9a7a0>
u_num u_den v_num v_den rel res
201 4 136 -133 (2, 0, -1, 2) 0
201 4 -1005 568 (5, 0, 1, 0) 0
201 4 4372152 935219 (-118, 82, 3, -14) 0
201 4 6210699 13897628 (-4, 201, 0, 0) 0
201 4 7919435 17426416 (91, 1, 11, -5) 0
-5 8 -1617 200 (8, 5, 0, 0) 0
-5 8 -477 692 (-3, 4, 3, 2) 0
-5 8 20824 2003 (8, 5, 0, 0) 0
-5 8 34272 -4885 (8, 5, 0, 0) 0
-5 8 36696 8687 (8, 5, 0, 0) 0
-5 8 398113 66200 (8, 5, 0, 0) 0
"""

"""
### How to use this to find new solutions

Once you run this on your 100 solutions, you will likely see one of two things:

1.  **The Constant Relation:** 
Many solutions might satisfy the *same* small relation (e.g., $2p - 3q + r = 0$). 
This means those solutions all lie on a specific **section** of the K3 surface. 

    **Strategy:** 
You can then substitute $r = 3q - 2p$ back into your quartic equation,
reducing the two-variable search $(u, v)$ into a one-variable search in $p/q$. 
This allows you to go to much higher heights.

2.  **The Bounded Height Relation:** 
If the relations are different but the coefficients $(m_1, ..., m_4)$ 
remain small, it suggests the solutions are generated by a specific 
**Automorphism** of the surface.

### Integrating with the Jacobian (The "Snapping" Method)
If the LLL residual is not zero but very small, you have a "near-miss." 
This is where you apply the **Newton-PSLQ** approach:

1.  **Numerical Refinement:** Take the $(u, v)$ that gave the small 
residual and use a Newton iteration on the quartic $D^2 = f(v)$ to 
find the "true" nearby rational point.
2.  **Jacobian Jump:** Once you have that "true" point (even if it's 
high height), map it to the Jacobian elliptic curve $E$.
3.  **The "Descent" Search:** Instead of searching for $u$ and $v$ 
integers, search for the **Rank** of $E$. If $E$ has rank >= 2$, 
use the generators to find the *other* points on that specific curve. 

This lattice-based view treats your list of 100 solutions as a 
**cloud of points** in a 4D space. LLL is essentially looking for the 
"planes" or "lines" that these points prefer to sit on. If you find 
a plane, your search space for the 101st solution shrinks by an 
entire dimension. 

**Which of your solutions have the largest $(u, v)$ heights?** 
If you provide a few more pairs, we can test if they share a 
common lattice relation.


I am using Piezas' method. But I am not familiar with the transformation 
$u to frac{u+1}{u-1}$ mentioned by Tito Piezas III

The transformation $u to frac{u+1}{u-1}$ is part of the automorphism 
group of the curves and surfaces you are studying. In the context of 
Tito Piezas III's work on Euler's quartic surface and the Fermat Quartic, 
this specific mapping is a linear fractional transformation 
(or Möbius transformation) of order 2.Here is how it functions and how 
you can use it to find the "hidden" smaller solutions:

1. The Nature of the Transformation
This transformation, let's call it $f(u)$, has a special property: 
it is its own inverse (an involution).
$$f(u) = frac{u+1}{u-1} implies f(f(u)) =
 frac{frac{u+1}{u-1} + 1}{frac{u+1}{u-1} - 1} =
   frac{(u+1) + (u-1)}{(u+1) - (u-1)} = 
   frac{2u}{2} = u$$
In the geometry of the surface, this map permutes the 48 lines. 
When you apply it to a $u$ value that generates an elliptic curve, 
it often transforms that curve into an isogenous one—meaning the 
new curve has the same "arithmetic soul" but different (often smaller) 
coefficients.

2. The "12 $u$ Values" ConnectionYou mentioned that every solution 
has 12 associated $u$ values. These values are generated by a group 
of transformations. The set typically includes:

Inversions: $u to 1/u$Negations: $u to -u$

The Piezas/Euler Shift: $u to frac{u+1}{u-1}$ 
and $u to frac{u-1}{u+1}$
By combining these (e.g., $u to -frac{u+1}{u-1}$), 
you get the full set of 12 rational points on the projective line 
$mathbb{P}^1$ that correspond to the same $A, B, C, D$ solution.

3. Why this helps you find "Huge" SolutionsThe "Huge" $(u, v)$ 
pairs you are finding are likely "huge" only because you are 
looking at them through one specific "window" (one $u$ value).

The Strategy:For your huge solution $u = 11846053/175812$:
Calculate $u_2 = frac{u+1}{u-1}$Calculate $u_3 = 1/u$
Calculate $u_4 = frac{1+u}{1-u}$...and so on for all 12 permutations.

The result:You will almost certainly find that one of these 12 
values is a very small fraction (like $5/8$ or $1/4$).
A "huge" $u$ near $1$ maps to a very large number via $frac{u+1}{u-1}$.
A "huge" $u$ near $1$ maps to a very small number via $u-1$.4. 
How to use it with Piezas' MethodPiezas often uses the parameter
 $k$ or $u$ to define the coefficients of the elliptic curve. 
 If your current search is struggling because the coefficients 
 $d_0, ..., d_4$ are becoming astronomical, switch the $u$.Instead of searching 
 $u in [10^6, 10^7]$, search for the transformed $u'$ in a smaller range. 
 If you find a solution at $u'$, you can map it back to the "huge" $u$ 
 using the inverse formula.Try this in your Sage code:Add a loop that tests 
 all 12 transformations of your current $u$ before giving up on it. 
 You might find that the "Gerbicz solutions" you think are out of reach 
 are actually sitting right next to the Elkies solutions, just viewed 
 through a different algebraic lens.Tip: In Piezas' notation, this 
 transformation is often linked to the "Double-Box" method for making 
 quartics into squares. If $u$ makes $f(u)$ a square, then 
 $frac{u+1}{u-1}$ will also likely satisfy a related condition on a dual curve.

-------------
Final note: Why I gave up:

I have set up a complex quadratic relation and put in strong
requirements for interation:
    # Quadratic Relation: 
    # Can still get trivial solutions when height of v much larger than u.
    M = Matrix(ZZ, [
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, p * 1 * K],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, p * p * K],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q * 1 * K],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q * q * K],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, p * q * K],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, r * 1 * K],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, r * r * K],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, s * 1 * K],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, s * s * K],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, r * s * K],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, p * r * K],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, p * s * K],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, q * r * K],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, q * s * K],
    ])
    
    # The first row of the LLL-reduced matrix is the shortest vector
    L = M.LLL()
    for i in range(L.nrows()):
        relation = L[0][:-1] # remove the 'K' column
        residual = L[0][-1]  # should be 0 if a perfect relation exists
        if 0 != residual: break # no good relation
        if not (2 < relation.hamming_weight() < 6):
            residual = -1
            break # too simple or too complex
        u_any = 0 < relation[:5].hamming_weight()
        v_any = 0 < relation[5:10].hamming_weight()
        inter = 0 < relation[10:].hamming_weight()
        if not (inter or (u_any and v_any)):
            residual = -1
            break # does not include both u and v
    return relation, residual

These are the results for all of the solutions in the table:
201 4 136 -133 -> (1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 1, 0, 0)
201 4 -1005 568 -> (1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
-5 8 -1617 200 -> (0, 1, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1)
-29 12 1865 132 -> (0, 0, 1, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
-29 12 -3333 107368 -> (-1, 4, 0, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
-9 20 -1425 412 -> (-1, -1, 1, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0)
-9 20 30080 6007 -> (1, -1, 0, 0, 0, 2, 0, -1, 0, 0, 0, 1, 0, 0)
-41 36 9360 2371 -> (-1, 0, 0, -3, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0)
-93 80 400 -37 -> (0, 1, 1, 0, -1, 0, 0, 0, -1, 1, 0, 0, 0, 0)
-93 80 -2433 920 -> (0, 0, -1, -2, 0, 0, 0, 1, 0, 0, 0, -1, 0, -1)
1865 132 30768 -57253 -> (5, 0, -2, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0)
136 -133 -1005 568 -> (2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 2, 0, 0)
400 -37 -2433 920 -> (0, -3, 0, 0, 1, 0, 0, 0, -1, 0, -1, 1, 0, 0)
12185 432 22529 2988 -> (0, 0, 0, 1, 4, 0, 0, 0, -1, 0, 0, 0, -1, -2)
1152 -2345 -15461 13160 -> (0, 0, 0, 0, 0, 3, -1, 0, 0, -1, -2, 0, 0, 0)

I don't see any repeats, so I don't know how to use the information.

Gemini AI offered some suggestions, but I don't think they are valid,
so give up on this method.
"""


def DEBUG(*args):
    set_trace()
    pass; pass; pass # opportumity to debug

if __name__ == "__main__":
    import sys
    command = eval(sys.argv[1])
    print(command)
    args = list(map(sage_eval, sys.argv[2:]))
    result = command(*args)
    print(result)
