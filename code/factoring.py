""" 
Select a fast partial factoring strategy
to be run on my MAC laptop with M1 chip.
The usage case is
    a^4 + b^4 = d^4 - c^4 for variables around 10^8
Define:
    t = a / 8, u = b / 8, v = (d - c) / 2, w = (d + c) / 2.
Then the equation becomes
    (t^4 y + u^4) * 2^9 = m = v * w * (v^2 + w^2) = v^3 * w + v * w^3
So, given (t, u), I want to fully factor unknown x = v or w.

I could use sympy or primefac packages if I wanted full factorization,
but I just want a high probability that one of the algebraic factors
v or w is completely factored.

When factoring the sum of 4th powers, the primes must be odd,
except for common factors of (t, u).
"""
from math import cbrt, gcd, isqrt, sqrt 
import random
from time import time
from logging_413 import V, IntFlag, parse_flags
from numpy import searchsorted
from probable_prime import is_prime, next_prime
from supplement.prime_sieve import PrimeProductTree
from utilities_413 import load_tree, factor_tree_gen

# Load and subset prime product trees
tree357_13 = load_tree(-2, 'tree357_13') # 3,5,7mod8, 14 levels, length 2**13, 3 thru 115_223
tree357_small = tree357_13[:7] # thru 421, sufficient for up to u = 100_000
tree1_17 = load_tree(-2, 'tree1_17') # 1mod8, 17 levels, length 2**16, 17 thru 3_684_641
tree1_small = tree1_17[:4] # thru 193, sufficient for up to u = 100_000

# Tune the size of the 1mod8 tree
#tree1 = tree1_17[:1] # just 17
#tree1 = tree1_17[:2] # plus 41
#tree1 = tree1_17[:3] # plus 73 thru 89
#tree1 = tree1_17[:4] # plus 97 thru 193
#tree1 = tree1_17[:5] # plus 233 thru 401
#tree1 = tree1_17[:6] # plus 409 thru 857
#tree1 = tree1_17[:7] # plus 881 thru 1801
#tree1 = tree1_17[:8] # plus 1873 thru 3889
#tree1 = tree1_17[:9] # plus 3929 thru 8377
#tree1 = tree1_17[:10] # plus 8513 thru 18_233
#tree1 = tree1_17[:11] # plus 18_257 thru 39_769
#tree1 = tree1_17[:12] # plus 39_841 thru 84_913
tree1 = tree1_17[:13] # plus 84_961 thru 181_777
#tree1 = tree1_17[:14] # plus 181_873 thru 389_057
#tree1 = tree1_17[:15] # plus 389_089 thru 823_777
#tree1 = tree1_17[:16] # plus 823_841 thru 1_747_489
#tree1 = tree1_17[:] # plus 1_747_489 thru 3_684_641

# Build list of all primes (357mod8 and 1mod8)
p1_inx = searchsorted(tree1_17[0], tree357_13[0][-1], 'right') # inx highest prime
primes_bounded = sorted([2] + tree357_13[0] + tree1_17[0][:p1_inx])[:len(tree357_13[0])]

# Tune the size of prime powers for p-1
#pp_len = 2**1
#pp_len = 2**2
#pp_len = 2**3
#pp_len = 2**4
#pp_len = 2**5
#pp_len = 2**6
pp_len = 2**7
#pp_len = 2**8
#pp_len = 2**9
#pp_len = 2**10
#pp_len = 2**11
#pp_len = 2**12
#pp_len = len(primes_bounded) # 2**13 = 8192

# Fill prime_powers with largest power q^k <= B
prime_powers = [1] * pp_len
B = primes_bounded[pp_len - 1] # max prime
for inx, q in enumerate(primes_bounded[:pp_len]):
    p_pow = q
    while p_pow <= B: p_pow *= q
    prime_powers[inx] = p_pow // q

# Chunk prime_powers into blocks for batched p-1 with backtracking
block_factor = 2**3
assert block_factor < pp_len
block_size = pp_len // block_factor
blocks = [1] * ((pp_len + block_size - 1) // block_size)
for pp_inx in range(pp_len):
    block_inx = pp_inx // block_size
    blocks[block_inx] *= prime_powers[pp_inx]

def p_minus_1_backtrack(n :int, vV: IntFlag=V.NONE) -> int | None:
    """
    Batched P-1 with backtracking when g == n.
 
    :param n (int): The odd integer to factor.
    :param vV (IntFlag): logging options
        
    :return (int or None): A non-trivial factor if found, otherwise None.    
    """
    a = 2 # standard base
   
    # Process blocks of primes to reduce GCD overhead
    for i, q in enumerate(blocks):
        last_a = a # Save checkpoint
        a = pow(a, q, n)
        g = gcd(a - 1, n)
        if 1 == g:
            continue
        if g < n:
            g = min(g, n // g)
            V.log(vV, V.F_P_1, f'p-1 block {i} found {g} in {n}')
            return g # Found a factor/fragment
        # OVERSHOOT: Multiple factors found. Backtrack.
        assert g == n
        a = last_a
        for pp_inx in range(i * block_size, pp_len):
            q = prime_powers[pp_inx]
            a = pow(a, q, n)
            g = gcd(a - 1, n)
            if 1 == g:
                continue
            if g < n:
                g = min(g, n // g)
                V.log(vV, V.F_P_1, f'p-1 backtrack pp_inx {pp_inx} found {g} in {n}')
                return g # Found a factor/fragment
            # OVERSHOOT again: no backtrack
            assert g == n
            V.log(vV, V.F_P_1, f'p-1 backtrack block {i} failed in {n}')
            return None
    V.log(vV, V.F_P_1, f'p-1 all blocks fail {n}')
    return None

# NOT USED EXCEPT FOR DEBUGGING
def p_minus_1(n :int, pp :int=prime_powers, vV: IntFlag=V.NONE) -> int | None:
    """Attempt to find a non-trivial factor of n using Pollard's p-1 algorithmn
 
    :param n (int): The odd integer to factor.
    :param vV (IntFlag): logging options
        
    :return (int or None): A non-trivial factor if found, otherwise None.    
    """
    a = 2
    for q in pp: # p^k       
        # Update a = a^(p^k) mod n
        a = pow(a, q, n)
        g = gcd(a - 1, n)
        if 1 == g:
            continue
        if g == n:
            V.log(vV, V.F_P_1, f'simple p-1 overshoot {n}')
            return None # Overshoot, no backtrack
        else:
            g = min(g, n // g)
            V.log(vV, V.F_P_1, f'simple p-1 found {g} in {n}')
            return g
    V.log(vV, V.F_P_1, f'simple p-1 fail {n}')
    return None


def pollard_rho_one(n :int, max_attempts :int=5, steps_1 :int=isqrt(2**17),
                    steps_mult :float=sqrt(2),  brent_m :int = 128,
                    vV: IntFlag=V.NONE) -> int | None:
    """
    Attempt to find ONE non-trivial factor of n using Pollard's rho algorithm. 
    
    :param n (int): The odd integer to factor.
    :param max_attempts (int): How many times to random restart.
    :param steps_1 (int): Initial number of steps based on expected factor
    :param steps_mult (float): step limit multiplier
        sqrt(2**17) * sqrt(2)**5 -> sqrt(2**21)
        362 * 6.656 -> ~1448
    :param brent_m (int): size of Brent loop
    :param vV (IntFlag): logging options
        
    :return (int or None): A non-trivial factor if found, otherwise None.
    """
    # Reduce initial steps_limit if n < expected, but not too far
    steps_limit = max(int(1.5 * brent_m), min(steps_1, 1 << (n.bit_length() // 4 + 2)))

    for attempt in range(max_attempts):
        if 0 < attempt:
            steps_limit = int(steps_limit * steps_mult)
        V.log(vV, V.F_RHO, f'rho attempt {attempt}, limit {steps_limit} on {n}')

        # Randomize starting position and polynomial constant
        y = random.randrange(2, n)
        c = random.randrange(1, n)
        g = 1
        r = 1
        steps = 0

        while g == 1 and steps < steps_limit:
            x = y
            # advance y by r steps
            for _ in range(r):
                y = (y * y + c) % n
                steps += 1

            k = 0
            while k < r and g == 1:
                ys = y
                q = 1
                bound = min(brent_m, r - k)
                for _ in range(bound):
                    y = (y * y + c) % n
                    # multiply differences into q modulo n
                    q = (q * (abs(x - y))) % n
                    steps += 1
                g = gcd(q, n)
                if g != 1:
                    break
                k += bound

            r <<= 1
            r = min(r, steps_limit - steps)

        if g == n: # failure in Brent's loop; try fallback linear search
            g = 1
            fallback_limit = max(steps_limit - steps, steps // 4)
            V.log(vV, V.F_RHO, f'rho fallback at {steps} limit {fallback_limit}')
            while g == 1 and steps < fallback_limit:
                ys = (ys * ys + c) % n
                g = gcd(abs(x - ys), n)
                steps += 1
 
        if g == 1 or g == n: # attempt failed
            continue

        assert 1 < g
        V.log(vV, V.F_RHO, f'rho found {g} in {n}; attempt {attempt} steps {steps}')
        return g # success

    return None # All attempts failed

def roots_and_rho(candidates: dict, factors: dict, others: dict, vV: IntFlag=V.NONE) -> None:
    """
    Partially factor cadidates into counts of prime factors and other factors.

    :param candidates (dict): unfactored candidates (may be 1, prime, or powers)
    :param factors (dict): prime factors to be updated
    :param others (dict): non-prime factors to be updated
    :param vV (IntFlag): logging options
    :return None
    """
    while 0 < len(candidates):
        V.log(vV, V.F_RHO, f'Pollard_rho factoring {candidates}')
        new_candidates = dict()
        for n, c in candidates.items():
            if 1 == n:
                continue
            # Check whether n is a perfect square or higher square
            while True:
                r = isqrt(n)
                if r * r == n:
                    n = r
                    c *= 2
                else: break
            # Check whether n is a perfect cube
            r = round(cbrt(n))
            if r**3 == n:
                n = r
                c *= 3
            # Check whether n is prime
            if is_prime(n):
                factors[n] = factors.get(n, 0) + c
                continue
            factor = pollard_rho_one(n, vV=vV)
            if factor is None:
                others[n] = others.get(n, 0) + c
                continue
            new_candidates[factor] = new_candidates.get(factor, 0) + c
            n //= factor
            new_candidates[n] = new_candidates.get(n, 0) + c
        candidates = new_candidates

def factor_common(t: int, u: int, vV: IntFlag=V.NONE) -> list:
    """
    Attempt to fully factor gcd(t, u)
    
    :param t:int smaller variable
    :param u:int larger variable
    :param vV:IntFlag logging options
    :return: list reduced t, u, factors
    """
    g = gcd(t, u)
    if 1 == g: return t, u, {}
    factors = dict() # counts of extracted prime factors
    t_r, u_r = t, u    
    if g % 2 == 0: # extract all 2s
        tz = (g & -g).bit_length() - 1
        factors[2] = factors.get(2, 0) + tz
        evens = 2**tz
        t_r //= evens
        u_r //= evens
        g //= evens
    for tree in (tree357_small, tree1_small): # extract small primes
        if 1 == g:
            break
        g2 = gcd(g, tree[-1][0])
        if is_prime(g2):
            c = 0
            while(g % g2 == 0):
                c += 1
                t_r //= g2
                u_r //= g2
                g //= g2
            assert 0 < c
            factors[g2] = factors.get(g2, 0) + c
        else:
            for p in factor_tree_gen(g2, tree):
                c = 0
                while(g % p == 0):
                    c += 1
                    t_r //= p
                    u_r //= p
                    g //= p
                assert 0 < c
                factors[p] = factors.get(p, 0) + c
    # NOTE: could check for integer power first, but none seen yet
    if 1 != g and is_prime(g):
        factors[g] = factors.get(g, 0) + 1
        t_r //= g
        u_r //= g
        g //= g
        assert 1 == g
    V.log(vV, V.F_COMMON, f'Common factors {factors}, reduced t, u {t_r}, {u_r}')
    assert g == 1, f'Unresolved common residue {g} in t_r {t_r}, u_r {u_r}' \
                f' reduced from t {t}, u {u}'
    return t_r, u_r, factors

def factor_tu(t: int, u: int, vV: IntFlag=V.NONE) -> list:
    """
    Partially factor t^4 + u^4 into prime factor counts and others.
   
    :param t:int smaller variable
    :param u:int larger variable
    :param vV:IntFlag logging options
    :return: list [common_factors, prime_factors, other_factors]
    """
    V.log(vV, V.FACTOR, f'Factor t, u: {t}, {u}')
    tf = time()
    t, u, common_factors = factor_common(t, u, vV)
    tc = time()
    factors = dict() # counts of extracted prime factors

    m = t**4 + u**4

    if m % 2 == 0: # extract all 2s
        tz = (m & -m).bit_length() - 1
        factors[2] = factors.get(2, 0) + tz
        evens = 2**tz
        m //= evens

    # extract p-1 primes
    candidates = []
    while True:
        if is_prime(m):
            factors[m] = factors.get(m, 0) + 1
            m = 1
            break
        #cand = p_minus_1(m, vV=vV)
        cand = p_minus_1_backtrack(m, vV=vV)
        if cand is None: break
        if is_prime(cand):
            factors[cand] = factors.get(cand, 0) + 1
        else: candidates.append(cand)
        m //= cand
    t1 = time()
    
    # extract 1mod8 primes
    if 1 < m: candidates.append(m)
    new_candidates = dict()
    for cand in candidates:
        V.log(vV, V.F_GCD, f'Feed {cand} to GCD tree')
        g = gcd(cand, tree1[-1][0])
        if 1 == g:
            new_candidates[cand] = new_candidates.get(cand, 0) + 1
            continue
        if is_prime(g):
            c = 0
            while(cand % g == 0):
                c += 1
                cand //= g
            assert 0 < c
            factors[g] = factors.get(g, 0) + c
        else:
            for p in factor_tree_gen(g, tree1):
                c = 0
                while(cand % p == 0):
                    c += 1
                    cand //= p
                assert 0 < c
                factors[p] = factors.get(p, 0) + c
        if 1 < cand:
            if is_prime(cand):
                factors[cand] = factors.get(cand, 0) + 1
            else:
                new_candidates[cand] = new_candidates.get(cand, 0) + 1
    tg = time()
    
    # Pollard Rho with checks whether components are roots or prime
    # Updates factors, others in place.
    others = dict()
    if 0 < len(new_candidates):
        roots_and_rho(new_candidates, factors, others, vV) 
    tr = time()
    V.log(vV, V.FACTOR, f'common {tc-tf}, p-1 {t1 -tc}, gcd {tg-t1}, rho {tr-tg}')

    V.log(vV, V.FACTOR, f'{t}^4 + {u}^4 -> {factors} and {others}')
    #import pdb; pdb.set_trace()
    return common_factors, factors, others


## Unit tests
# ----------------------------------------------------------------------------
import unittest

class TestFactoringSearch412Fast(unittest.TestCase):

    def test_dummy_fast(self):
        self.assertTrue(True)

### END unittests 

## --- Main Section ---
import argparse
import sys

def main(argv=None):
    """Command-line dispatcher"""
    parser = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers(dest='command')

    # p_1 command
    p_p_1 = subparsers.add_parser('p_1', help='Factor with p_minus_one')
    p_p_1.add_argument('n', type=str, help='Python expression for number to factor')
    p_p_1.add_argument('-v', '--vV', type=parse_flags, default=V.NONE,
        help="Set verbosity levels (e.g., 'OUTER', 'OUTER|PROCESS', 'DEBUG')")

    # rho1 command
    p_rho1 = subparsers.add_parser('rho1', help='Factor with pollard_rho_one')
    p_rho1.add_argument('n', type=str, help='Python expression for number to factor')
    p_rho1.add_argument('max_attempts', type=str, nargs='?', default='5',
                        help='Python expression for max_attempts')    
    p_rho1.add_argument('-v', '--vV', type=parse_flags, default=V.NONE,
        help="Set verbosity levels (e.g., 'OUTER', 'OUTER|PROCESS', 'DEBUG')")

    # factor_tu command
    p_factor_tu = subparsers.add_parser('factor_tu', help='Factor with factor_tu')
    p_factor_tu.add_argument('t', type=str, help='Python expression for t')
    p_factor_tu.add_argument('u', type=str, help='Python expression for u')
    p_factor_tu.add_argument('-v', '--vV', type=parse_flags, default=V.NONE,
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

    if args.command == 'p_1':
        n = eval(args.n)
        start = time()
        g = p_minus_1_backtack(n, vV=args.vV)
        elapsed = time() - start
        print(f'Found {g}, Elapsed: {elapsed:.6f}s')
        return

    if args.command == 'rho1':
        n = eval(args.n)
        max_attempts = eval(args.max_attempts)
        start = time()
        g = pollard_rho_one(n, max_attempts=max_attempts, vV=args.vV)
        elapsed = time() - start
        print(f'Found {g}, Elapsed: {elapsed:.6f}s')
        return

    if args.command == 'factor_tu':
        import pdb; pdb.set_trace()
        t = eval(args.t)
        u = eval(args.u)
        start = time()
        results = factor_tu(t, u, vV=args.vV)
        elapsed = time() - start
        print(f'Found {results}, Elapsed: {elapsed:.6f}s')
        return

    # No command, so run unit tests
    return unittest.main(argv=[sys.argv[0]] + rest)

if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
