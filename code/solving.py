"""
As explained in README_code.md, we are given the factors of
m in the equations
    t^4 + u^4 = m; 2^9 * m = x * (y^3 + x^2 * y)
We select small subsets of the identified factors with product x .
If we can solve the monotonically increasing cubic function
    f(y; x) = y^3 + x^2 * y - 2^9 * m / x
for the root y, then x and y solve the equation.
"""
from math import cbrt, prod
from time import time
from logging_413 import V, IntFlag, parse_flags

def build_partition_bank(common: dict, factors: dict, others: dict) -> tuple:
    """Build a list of all factors from the three dicts for partitioning.
    
    All of the 2s are assigned to the left_product.
    Any 2s in factors are multiplied by the default 2**9 in the equation.
        This is not actually necessary, but simplifies the search space.
        Furthermore, the x or y containing the 2s is more likely to be fully factored.
    All of the `common` factors are raised to 4th power for distribution.
    The remaining factors and their powers are partitioned in all possible ways.
    The `others` factors are chosen for the left product last because they are
        more likely to be in the (xy^3 + yx^3) term.
    
    :param common: dict of common factors of t and u
    :param factors: dict of prime factors of m
    :param others: dict of other factors of m
 
    :return (list): [total_product (including 2**9), 
                     twos_product, 
                     [list of fragments for partitioning] ]

    Example:
        common = {5: 1); factors = {2: 10, 26777: 1, 102481: 1}; others = {4216393: 1}
        -> [7231466424844150625, 2**10, [5**4, 26777, 102481, 4216393] ]
    """
    # Handle 2s and initial left/right products
    twos_product = 2**9
    if 2 in common:
        twos_product *= 2**(4 * common[2])
        del common[2]
    if 2 in factors:
        twos_product *= 2**factors[2]
        del factors[2]

    # Fill bank of factors to partition
    bank_len = sum(common.values()) + sum(factors.values()) + sum(others.values())
    bank = [1]*bank_len
    bank_inx = 0
    for p, count in common.items():
        for _ in range(count):
            bank[bank_inx] = p**4
            bank_inx += 1
    for d in (factors, others):
        for p, count in d.items():
            for _ in range(count):
                bank[bank_inx] = p
                bank_inx += 1
    assert bank_inx == bank_len

    m512 = prod(bank) * twos_product # total product of all factors including 2**9
    return m512, twos_product, bank

def partition_fragments(m512: int, twos_product: int, bank: list,
                        range=range, len=len) -> int:
    """Yield left_product where left_product < right_product
    and the products come from partitioning `factors` into two groups and multiplying.
    
    Two twos_product is always assigned to the left_product.
        This is not actually necessary, but simplifies the search space.
        Furthermore, the x or y containing the 2s is more likely to be fully factored.
    The fragments in bank are partitioned in all possible ways.
    
    :param m512: total product of all factors including 2**9   
    :param twos_product: product of all 2s
    :param bank: list of fragments for partitioning 
 
    ;yield: left_product
    
    Notes:
    - The number of partitions is O(2^n) for `n = len(factors)`; use with
      care for large `n`.

    Example:
        7231466424844150625, 2**10, [5**4, 26777, 102481, 4216393]
        -> 2**10,
           2**10 * 26777,
           2**10 * 5**4 * 26777,
           2**10 * 102481,
           2**10 * 5**4 * 102481},
           ....
    """
    init_left = twos_product
    init_right = m512 // twos_product
    if init_left **2 >= init_right: return None

    # Generate partitions
    left_product_set = set()
    add_to_set = left_product_set.add # localize
    partition_index = -1 # Each bit flags whether corresponding factor is in left_product
    max_index = 2 ** len(bank) - 1 # all factor inclusion flags on
    while partition_index < max_index:
        # Bump partition index and use to select factors
        partition_index += 1
        left_product = init_left
        right_product = init_right
        for i in range(len(bank)):
            if partition_index & 2**i:
                v = bank[i]
                left_product *= v
                right_product //= v
        # Skip if reached max left product
        if left_product **2 >= right_product:
            continue
        # Avoid duplicates
        if left_product in left_product_set:
            continue
        add_to_set(left_product)
        # Yield partition
        yield left_product

    return None

def solve_cubic_by_cardano(x: int, c: int, cbrt=cbrt, round=round) -> int | None:
    """Attempt to find an integer root of 
        y^3 + x^2 * y - c = 0
    for y using Cardono's formula for the depressed cubic.
        y^3 + p * y + q = 0
    so p = x^2, q = -c.

    Since x^2 > 0, the cubic is strictly increasing (monotonic), and hence
    there is at most one real root.
    Args:
        x (int): sqrt of linear coefficient.
        c (int): constant term.
    
    Returns:
        int | None: The integer y or None.
    """
    p = x * x
    q = -c
    discriminant = (q / 2) ** 2 + (p / 3) ** 3
    if discriminant < 0:
        return None  # No real roots

    sqrt_discriminant = discriminant ** 0.5
    u_cubed = -q / 2 + sqrt_discriminant
    v_cubed = -q / 2 - sqrt_discriminant

    u = round(cbrt(u_cubed)) if u_cubed >= 0 else -round(cbrt(-u_cubed))
    v = round(cbrt(v_cubed)) if v_cubed >= 0 else -round(cbrt(-v_cubed))

    y = u + v
    if y < 0:
        return None

    # Verify that y is an integer root
    if y**3 + p * y == c:
        return y
    return None

def solve_xm(x: int, m: int, vV=V.NONE) -> list|None:
    """Given x and m in the equation
        t^4 + u^4 = m; 2^9 * m = x * (y^3 + x^2 * y)
    search for integer solution x, y and return sorted pair v, w.

    :param x: int value for x. Must be divisible by 2**9
    :param m: int value for m. Must be divisible by x // 2**9
    :param vV:IntFlag logging options
    :return: solution [v, w], or None if none found

    EXAMPLE (Solution #0)
    For t, u = 11975, 51820: m = 2^-9 * v * w * (v^2 + w^2)
        m = 102481 * (320000 = 2**9 * 5**4) * 112902355361 // 2**9
        m = 7231466424844150625
        given x = 320000 or 102481
        -> [v, w] = [102481, 320000]

    EXAMPLE (Solution #1)
    For t, u = 335305, 2349595: m = 2^-9 * v * w * (v^2 + w^2)
        m = 2625017 * (17990656 = 2**10 * 17569) * 330554417560625 // 2**9
        m = 30489627906502870450351250
        given x = 2625017 or 17990656
        -> [v, w] = [2625017, 17990656]
    """
    assert (2**9 * m) % x == 0, "2^9 * m must be divisible by x"
    c =  2**9 * m // x
    return solve_cubic_by_cardano(x, c)

def solve_factors(common: dict, factors: dict, others: dict, vV: IntFlag=V.NONE,
                  partition_fragments=partition_fragments, 
                  solve_cubic_by_cardano=solve_cubic_by_cardano) -> list|None:
    """Given factors of m in the equations
        t^4 + u^4 = m; 2^9 * m = x * (y^3 + x^2 * y)
    search for integer solution x, y and return sorted pair v, w.

    :param common: dict of common factors of t and u
    :param factors: dict of prime factors of m
    :param others: dict of other factors of m
    :param vV:IntFlag logging options
    :return: solution [v, w], or None if none found

    EXAMPLE (Solution #0)
    Factors for t, u = 11975, 51820:
        common {5: 1}
        factors {26777: 1, 102481: 1, 4216393: 1}
        others {}
        -> [v, w] = [102481, 320000]

    EXAMPLE (Solution #1)
    Factors for t, u = 335305, 2349595:
        common {5: 1}
        factors {2: 1, 97: 1, 353: 1, 761: 1, 17569: 1,20297: 1, 2625017: 1}
        others {}
        -> [v, w] = [2625017, 17990656]
    """
    m512, twos_product, bank = build_partition_bank(common, factors, others)
    for i, x in enumerate(partition_fragments(m512, twos_product, bank)):
        y = solve_cubic_by_cardano(x, m512 // x)
        if y is None:
            continue
        v, w = (x, y) if x < y else (y, x)
        V.log(vV, V.SOLVE, f"Found solution v={v}, w={w} in partition {i}")
        return [v, w]
    return None

def solve_factors_all(common: dict, factors: dict, others: dict, vV: IntFlag=V.NONE,
                  partition_fragments=partition_fragments, 
                  solve_cubic_by_cardano=solve_cubic_by_cardano) -> list|None:
    """Same as solve_factors, but does not stop if a factor is found.
    """
    solution = None
    m512, twos_product, bank = build_partition_bank(common, factors, others)
    for i, x in enumerate(partition_fragments(m512, twos_product, bank)):
        y = solve_cubic_by_cardano(x, m512 // x)
        if y is None:
            continue
        v, w = (x, y) if x < y else (y, x)
        V.log(vV, V.SOLVE, f"Found solution v={v}, w={w} in partition {i}")
        assert solution is None # expect only one solution for given factors
        solution = [v, w]
    return solution



## Unit tests
# ----------------------------------------------------------------------------
import unittest

class TestSolveXMFast(unittest.TestCase):
    def test_solution0(self):
        r = solve_xm(320_000, 7231466424844150625)
        self.assertEqual(r, 102_481)
        r = solve_xm(102_481, 7231466424844150625)
        self.assertEqual(r, 320_000)
        r = solve_xm(320_000, 7231466424844150625 * 2)
        self.assertIsNone(r)
        r = solve_xm(320_000, 7231466424844150625 * 3)
        self.assertIsNone(r)

### END unittests 

## --- Main Section ---
import argparse
import sys
import timeit

def main(argv=None):
    """Command-line dispatcher"""
    parser = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers(dest='command')

    # solve_tu command
    p_solve_tu = subparsers.add_parser('solve_tu', help='Solve with solve_tu')
    p_solve_tu.add_argument('t', type=str, help='Python expression for t')
    p_solve_tu.add_argument('u', type=str, help='Python expression for u')
    p_solve_tu.add_argument('-a', '--all', action='store_true', help="Report all")
    p_solve_tu.add_argument('-v', '--vV', type=parse_flags, default=V.NONE,
        help="Set verbosity levels (e.g., 'OUTER', 'OUTER|PROCESS', 'DEBUG')")

    # solve_xm command
    p_solve_xm = subparsers.add_parser('solve_xm', help='Solve with solve_xm')
    p_solve_xm.add_argument('x', type=str, help='Python expression for x')
    p_solve_xm.add_argument('m', type=str, help='Python expression for m')
    p_solve_xm.add_argument('-v', '--vV', type=parse_flags, default=V.NONE,
        help="Set verbosity levels (e.g., 'OUTER', 'OUTER|PROCESS', 'DEBUG')")

    """ parse known args and leave the rest for unittest
        -v	--verbose	Show individual test names and results.
        -q	--quiet	Show minimal output.
        -f	--failfast	Stop on first failure.
        -b	--buffer	Hide stdout/stderr for passing tests.
        -c	--catch	Graceful interrupt with Ctrl-C.
        -k	(None)	Run tests matching a pattern
    """
    args, rest = parser.parse_known_args(argv)
    if args.command == 'solve_tu':
        from factoring import factor_tu, return_saved_times
        t = eval(args.t)
        u = eval(args.u)
        common, factors, others = factor_tu(t, u, vV=args.vV)
        if args.all:
            results = solve_factors_all(common, factors, others, vV=args.vV)
            elapsed = timeit.timeit(lambda: 
                solve_factors_all(common, factors, others), number=1000) / 1000
        else:
            results = solve_factors(common, factors, others, vV=args.vV)
            elapsed = timeit.timeit(lambda: 
                solve_factors(common, factors, others), number=1000) / 1000
        print(f'timeit: {elapsed:.8f}s')
        print(f'Found {results}')
        if args.vV and args.vV & V.F_DIAG:
            print(return_saved_times())
        return
 
    if args.command == 'solve_xm':  
        x = eval(args.x)
        m = eval(args.m)
        start = time()
        results = solve_xm(x, m, vV=args.vV)
        elapsed = time() - start
        print(f'Found {results}, Elapsed: {elapsed:.6f}s')
        return
    # No command, so run unit tests
    return unittest.main(argv=[sys.argv[0]] + rest) 

if __name__ == "__main__":
    main()
