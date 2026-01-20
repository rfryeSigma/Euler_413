""" 
Select a fast partial factoring strategy
to be run on my MAC laptop with M1 chip.
The usage case is
    a^4 + b^4 = d^4 - c^4 for variables around 10^8
Define:
    t = a / 8, u = b / 8, v = (d - c) / 2, w = (d + c) / 2.
Then the equation becomes
    (t^4 y + u^4) * 2^9 = m = v * w * (v^2 + w^2) = v^3 * w + v * w^3
So, given (t, u), Try to fully factor at least on unknown x = v or w.

I could use sympy or primefac packages if I wanted full factorization,
but I just want a high probability that one of the algebraic factors
v or w is completely factored.

When factoring the sum of 4th powers, the primes must be 2s or 1mod8,
except for common factors of (t, u).
"""
from math import cbrt, exp, isqrt
from gmpy2 import gcd, is_square, is_prime
import random
from time import time
from logging_413 import V, IntFlag, parse_flags
from numpy import searchsorted
from utilities_413 import load_tree, factor_tree_gen

# Load and subset prime product trees
tree357_13 = load_tree(-2, 'tree357_13') # 3,5,7mod8, 14 levels, length 2**13, 3 thru 115_223
tree357_small = tree357_13[:7] # thru 421, sufficient for up to u = 100_000

level = 12 # select level (counting from primes at level 0) for processing tree1
tree1 = [[]]* (level + 1)
for i in range(level + 1): # truncate the top levels of prime factor tree
    tree1[i] = load_tree(i, 'tree1_18') # 1mod8, 18 levels, length 2**17, 17 thru 7_766_713
tree1_small = tree1[:4] # thru 193, sufficient for up to u = 100_000

######## TUNE THIS FOR OPTIMAL FACTORS, AVG TIME #########
num_indexes = 20 # number of indexes at selected level to scan for gcd test and factor
assert num_indexes <= len(tree1[-1]) # 32

# Build list of primes (both 357mod8 and 1mod8)
assert tree1[0][-1] > tree357_13[0][-1]
# Find index in 1mod8 primes of highest prime in 357modd8 primes
p1_inx = searchsorted(tree1[0], tree357_13[0][-1], 'right') # inx highest prime
# Merge the primes
primes_bounded = sorted([2] + tree357_13[0] + tree1[0][:p1_inx]) # len 10_891; ... 115_201, 115_211, 115_223]
# Tune the size of prime powers for p-1
#pp_len = 2**1 # thru 3
#pp_len = 2**2 # thru 7
#pp_len = 2**3 # thru 19
#pp_len = 2**4 # thru 53
#pp_len = 2**5 # thru 131
#pp_len = 2**6 # thru 311
#pp_len = 2**7 # thru 719
#pp_len = 2**8 # thru 1619
pp_len = 2**9 # thru 3671
#pp_len = 2**10 # thru 8161

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

def p_minus_1_backtrack(n :int,
                        pow=pow, gcd=gcd) -> int | None:
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
            #V.log(vV, V.F_P_1, f'p-1 block {i} found {g} in {n}')
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
                #V.log(vV, V.F_P_1, f'p-1 backtrack pp_inx {pp_inx} found {g} in {n}')
                return g # Found a factor/fragment
            # OVERSHOOT again: no backtrack
            assert g == n
            #V.log(vV, V.F_P_1, f'p-1 backtrack block {i} failed in {n}')
            return None
    #V.log(vV, V.F_P_1, f'p-1 all blocks fail {n}')
    return None

# NOT USED EXCEPT FOR DEBUGGING
def p_minus_1(n :int, pp :list=prime_powers,
              pow=pow, gcd=gcd, min=min) -> int | None:
    """Attempt to find a non-trivial factor of n using Pollard's P-1 algorithmn
 
    :param n (int): The odd integer to factor.
    :param pp (list): list of primes raised to powers up to limit.
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
            #V.log(vV, V.F_P_1, f'simple p-1 overshoot {n}')
            return None # Overshoot, no backtrack
        else:
            g = min(g, n // g)
            #V.log(vV, V.F_P_1, f'simple p-1 found {g} in {n}')
            return g
    #V.log(vV, V.F_P_1, f'simple p-1 fail {n}')
    return None

# Globals for pollard rho routines
max_attempts = 3
brent_m = 2**8
steps_mult = 1.6
def pollard_rho_probability(p: int, max_attempts: int=max_attempts + 1,
        steps_limit: int=3*brent_m, steps_mult: int=steps_mult) -> None:
    """Calculate the probability the pollard_rho_one will fail
    to find prime p.
    """
    prob_fail_cumulative = 1.0
    for attempt in range(max_attempts):
        steps_limit = int(steps_limit * steps_mult ** attempt)
        prob_fail = exp(- steps_limit * steps_limit / (2.0 * p))
        prob_fail_cumulative *= prob_fail
        print(f'attempt {attempt}, steps {steps_limit}'
              f', prob_fail {prob_fail:.2f}'
              f', prob_fail_cumulative {prob_fail_cumulative}:.2f')

def pollard_rho_one(n: int, max_attempts: int=max_attempts, brent_m: int=brent_m,
        steps_1: int=3*brent_m, steps_mult: int=steps_mult, vV: IntFlag=V.NONE,
        max=max, min=min, log=V.log, gcd=gcd) -> int | None:
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
    steps_limit = max(3 * brent_m, min(steps_1, n.bit_length() * 100))

    for attempt in range(max_attempts):
        steps_limit = int(steps_limit * steps_mult ** attempt)
        log(vV, V.F_RHO, f'rho attempt {attempt}, limit {steps_limit} on {n}')

        # Randomize starting position and polynomial constant
        y = random.randrange(2, n)
        c = random.randrange(1, n)
        g = 1
        r = 1
        steps = 0

        while g == 1 and steps < steps_limit:
            x = y
            k = 0
            while k < r and g == 1:
                checkpoint_y = y 
                q = 1
                bound = min(brent_m, r - k)
                for _ in range(bound):
                    y = (y * y + c) % n
                    q = (q * (x - y)) % n # Can skip abs here; gcd will fix
                    steps += 1
                g = gcd(q, n)
                if g != 1:
                    break
                k += bound

            if g == n: # failure in Brent's loop; try fallback at checkpoint        
                g = 1
                y_backtrack = checkpoint_y
                while g == 1:
                    y_backtrack = (y_backtrack * y_backtrack + c) % n
                    g = gcd(x - y_backtrack, n) # can skip abs here; gcd will fix
                if g == n:
                    break # Rare backtrack failure. Try again.
            r <<= 1
            r = min(r, steps_limit - steps)

        if g == 1 or g == n: # attempt failed
            continue

        assert 1 < g
        log(vV, V.F_RHO, f'rho found {g} in {n}; attempt {attempt} steps {steps}')
        return g # success

    return None # All attempts failed

def roots_and_rho(candidates: dict, factors: dict, others: dict, vV: IntFlag=V.NONE,
                  len=len, log=V.log, isqrt=isqrt, is_prime=is_prime, round=round,
                  cbrt=cbrt, pollard_rho_one=pollard_rho_one) -> None:
    """
    Partially factor cadidates into counts of prime factors and other factors.

    :param candidates (dict): unfactored candidates (may be 1, prime, or powers)
    :param factors (dict): prime factors to be updated
    :param others (dict): non-prime factors to be updated
    :param vV (IntFlag): logging options
    :return None
    """
    factors_get = factors.get
    others_get = others.get
    while 0 < len(candidates):
        log(vV, V.F_RHO, f'Pollard_rho factoring {candidates}')
        candidates_items = candidates.items
        candidates_get = candidates.get
        new_candidates = dict()
        new_candidates_get = new_candidates.get
        for n, c in candidates_items():
            if 1 == n:
                continue
            # Check whether n is a perfect square or cube
            if is_square(n):
                n = isqrt(n)
                c *= 2
            if round(cbrt(n))**3 == n:
                n= round(cbrt(n))
                c *= 3
            assert not is_square(n) and round(cbrt(n))**3 != n
            # Check whether n is prime
            if is_prime(n):
                factors[n] = factors_get(n, 0) + c
                continue
            if vV&V.RHO1: # call pollard_rho_one
                factor = pollard_rho_one(n, vV=vV)
                if factor is None:
                    others[n] = others_get(n, 0) + c
                    continue
                new_candidates[factor] = new_candidates_get(factor, 0) + c
                n //= factor
                new_candidates[n] = new_candidates_get(n, 0) + c
            else:
                others[n] = others_get(n, 0) + c
                continue
        candidates = new_candidates
 
def factor_common(t: int, u: int,
                  gcd=gcd, is_prime=is_prime, factor_tree_gen=factor_tree_gen,
                  log=V.log) -> tuple:
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
    factors_get = factors.get #localize
    t_r, u_r = t, u    
    if g % 2 == 0: # extract all 2s
        tz = (g & -g).bit_length() - 1
        factors[2] = factors_get(2, 0) + tz
        t_r >>= tz
        u_r >>= tz
        g >>= tz
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
            factors[g2] = factors_get(g2, 0) + c
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
    #log(vV, V.F_COMMON, f'Common factors {factors}, reduced t, u {t_r}, {u_r}')
    assert g == 1, f'Unresolved common residue {g} in t_r {t_r}, u_r {u_r}' \
                f' reduced from t {t}, u {u}'
    return t_r, u_r, factors

def factor_tu(t: int, u: int, vV: IntFlag=V.NONE,
              log=V.log, factor_common=factor_common, time=time,
              is_prime=is_prime, gcd=gcd, len=len,
              p_minus_1_backtrack=p_minus_1_backtrack,
              factor_tree_gen=factor_tree_gen,
              ) -> tuple:
    """
    Partially factor t^4 + u^4 into prime factor counts and others.
   
    :param t:int smaller variable
    :param u:int larger variable
    :param vV:IntFlag logging options
    :return: list [common_factors, prime_factors, other_factors]
    """
    log(vV, V.FACTOR, f'Factor t, u: {int(t):_}, {int(u):_}')
    tf = time()
    t, u, common_factors = factor_common(t, u)
    tc = time()
    factors = dict() # counts of extracted prime factors
    factors_get = factors.get #localize

    m = t**4 + u**4
    n = m # remainder of m after extracting factors

    if n % 2 == 0: # extract all 2s
        tz = (n & -n).bit_length() - 1
        factors[2] = factors_get(2, 0) + tz
        n >>= tz
    else: tz = 0

    # extract p-1 primes
    candidates = []
    while n > 1:
        if is_prime(n):
            factors[n] = factors_get(n, 0) + 1
            n = 1
            break
        f = p_minus_1_backtrack(n)
        if f is None: break
        if is_prime(f):
            factors[f] = factors.get(f, 0) + 1
        else: candidates.append(f)
        n //= f
    t1 = time()

    # Extract 1mod8 primes
    if 1 < n: candidates.append(n)
    p1_frags = len(candidates) + sum(factors.values()) - tz
    for index, test_gcd in enumerate(tree1[-1][:num_indexes]):
        new_candidates = dict()
        for n in candidates:
            if is_prime(n):
                factors[n] = factors_get(n, 0) + 1
                continue
            log(vV, V.F_GCD, f'Feed candidate {n} to {index} GCD tree')
            g = gcd(n, test_gcd)
            if 1 == g:
                new_candidates[n] = new_candidates.get(n, 0) + 1
                continue
            if is_prime(g):
                pc = 0
                while(n % g == 0):
                    pc += 1
                    n //= g
                factors[g] = factors_get(g, 0) + pc
            else:
                for p in factor_tree_gen(g, tree1, level, index, _func=factor_tree_gen):
                    pc = 0
                    while (n % p == 0):
                        pc += 1
                        n //= p
                    factors[p] = factors.get(p, 0) + pc
            if 1 < n:
                if is_prime(n):
                    factors[n] = factors.get(n, 0) + 1
                else:
                    new_candidates[n] = new_candidates.get(n, 0) + 1
        if 0 == len(new_candidates): break
        candidates = new_candidates
    tg = time()

    # Pollard Rho with checks whether components are roots or prime
    # Updates factors, others in place.
    others = dict()
    if 0 < len(new_candidates):
        roots_and_rho(new_candidates, factors, others, vV) 
    tr = time()

    # Verify factorization
    n = 1
    for p, c in factors.items():
        n *= p**c
    for p, c in others.items():
        n *= p**c
    assert n == m
    
    if vV and vV & V.F_DIAG:
        save_times(tf, tc, t1, tg, tr, p1_frags, factors, others)

    log(vV, V.FACTOR, f'{t}^4 + {u}^4 -> {factors} and {others}')
    return common_factors, factors, others

saved_common  = 0.0 # Elapsed time on common factors
saved_p1      = 0.0 # Elapsed time on p-1 method
saved_gcd     = 0.0 # Elapsed time on factor tree gcd method
saved_rho     = 0.0 # Elapsed time on rho method
saved_count   = 0 # Number of saved times
saved_frags = 0 # Number of fragments found by p-1
saved_factors = 0 # Number of factors
saved_others  = 0 # NUmber of others
def save_times(tf: float, tc: float, t1: float, tg: float, tr: float,
               p1_frags: dict, factors: dict, others: dict) -> None:
    """Save cumulative factor count and elapsed time on components"""
    global saved_common, saved_p1, saved_gcd, saved_rho, saved_count, \
        saved_frags, saved_factors, saved_others
    saved_common  += tc-tf
    saved_p1      += t1-tc
    saved_gcd     += tg-t1
    saved_rho     += tr-tg
    saved_count   += 1
    saved_frags  += p1_frags
    saved_factors += sum(factors.values())
    saved_others  += sum(others.values())

def return_saved_times() -> None:
    """Report averaged factor component times"""
    if saved_count == 0:
        print('No saved times')
        return
    total = (saved_common + saved_p1 + saved_gcd + saved_rho)
    total_100 = total / 100.0
    c = float(saved_count)
    return(f'count {saved_count:_}:' \
          f'\np1_frags {saved_frags / c:.6g}, factors {saved_factors / c:.6g},' \
          f' others {saved_others / c:.6g}, avg time {total / c:.3g}' \
          f'\ncommon {saved_common / c:.3g} ({saved_common / total_100:.02f}%)' \
          f', p-1 {saved_p1 / c:.3g} ({saved_p1 / total_100:.02f}%)' \
          f', gcd {saved_gcd / c:.3g} ({saved_gcd / total_100:.02f}%)' \
          f', rho {saved_rho / c:.3g} ({saved_rho / total_100:.02f}%)')

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
import timeit

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
    p_factor_tu.add_argument('-e', '--elapsed', action='store_true', 
                             help="Calculate elapsed with timeit")
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
        g = p_minus_1_backtrack(n, vV=args.vV)
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
        t = eval(args.t)
        u = eval(args.u)
        results = factor_tu(t, u, vV=args.vV)
        if args.elapsed:
            elapsed = timeit.timeit(lambda: factor_tu(t, u), number=100) / 100
            print(f'timeit: {elapsed:.8f}s')
        print(f'Found {results}')
        if args.vV and args.vV & V.F_DIAG:
            print(return_saved_times())
        return

    # No command, so run unit tests
    return unittest.main(argv=[sys.argv[0]] + rest)

if __name__ == "__main__":
    main()
