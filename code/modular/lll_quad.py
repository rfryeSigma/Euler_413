"""
Both of thsese suggestions are silly. 
They underestimate the chances of finding squares.

The suggestion to use LLL on the form (x^2)^2 + (y^2)^2 + (z^2)^2 = 1 
is a clever "linearization" tactic. It treats the problem not as a 
quartic equation, but as a quadratic form in the squares of the variables.

Let X = x^2, Y = y^2, Z = z^2. The equation becomes:
X^2 + Y^2 + Z^2 = (w^2)^2
N is a sum of three squares if and only if it is not of the form 4^k(8m + 7). 
There are far more integer points on a sphere than on a quartic surface.
"""
def lll_quadric_filter(limit_W, precision=40):
    # This is a conceptual sketch of the 'Near Miss' lattice
    # We want X^2 + Y^2 + Z^2 - W^2 = 0
    # Look for points near the sphere
    K = 10**precision
    for _ in range(100):
        # Pick a random direction on the sphere
        target = vector(RR, [random(), random(), random()])
        target = target / target.norm()
        
        # Build LLL matrix to find integers close to this direction
        M = Matrix(ZZ, [
            [1, 0, 0, floor(K * target[0])],
            [0, 1, 0, floor(K * target[1])],
            [0, 0, 1, floor(K * target[2])],
            [0, 0, 0, K]
        ])
        reduced = M.LLL()
        for row in reduced:
            X, Y, Z = abs(row[0]), abs(row[1]), abs(row[2])
            # Check if they are perfect squares
            if is_square(X) and is_square(Y) and is_square(Z):
                # Check the Fermat condition
                if X**2 + Y**2 + Z**2 == is_square(X**2 + Y**2 + Z**2):
                    return (sqrt(X), sqrt(Y), sqrt(Z))
    return None

"""
Combine this with 100 known solutions.
Square (a, b, c, d) to get (X, Y, Z, W).
Check if these squared points lie on a sub-lattice of the sphere 
    (do they satisfy mX + nY + kZ = 0 )
If they do, use that sub-lattice in the LLL search. 
This is significantly faster than searching the whole sphere.
"""
