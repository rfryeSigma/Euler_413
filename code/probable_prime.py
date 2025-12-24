# Miller-Rabin Probabilistic Primality Test

import random
import os


## The Miller-Rabin Primality Test
# ----------------------------------------------------------------------------
DETERMINISTIC_BASES_64 = [2, 325, 9375, 28178, 450775, 9780504, 1795265022]

def power(a, d, n):
    """Fast modular exponentiation wrapper using Python's built-in pow.

    Kept as a small wrapper so existing callers can continue to call
    `power(a, d, n)` while benefiting from the C-optimized `pow`.
    """
    return pow(a, d, n)

def is_prime_miller_rabin(n, k=8, deterministic: bool = False):
    """
    Performs the Miller-Rabin primality test.

    Args:
        n (int): The number to test for primality.
        k (int): The number of bases (iterations) to test. Higher k means
                 higher confidence. Default k=8 is generally sufficient for
                 ~10-digit numbers while keeping the test fast. Use larger
                 `k` for bigger inputs or higher confidence.
        deterministic (bool): If True and `n` is within 64-bit range, use a
                 deterministic set of bases that makes Miller-Rabin
                 deterministic for 64-bit integers. Otherwise bases are
                 chosen pseudo-randomly.

    Returns:
        bool: True if n is probably prime, False if n is definitely composite.
    """
    
    # 1. Handle edge cases
    if n <= 1:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False

    # 2. Factor n-1 into (2^s) * d, where d is odd
    # n-1 = 2^s * d
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    
    # 3. Choose bases. Use deterministic base set for 64-bit numbers when
    # requested; otherwise pick `k` pseudo-random bases. Precompute bases to
    # avoid calling the RNG inside the test loop.
    bases = []
    # Use deterministic bases for 64-bit integers when requested
    if deterministic and n < (1 << 64):
        # Use the known deterministic set; if k < len(bases) then truncate.
        bases = DETERMINISTIC_BASES_64[:k] if k > 0 else DETERMINISTIC_BASES_64[:]
    else:
        if k <= 0:
            bases = []
        else:
            # If the candidate space is larger than k, sample without replacement
            # from the integer interval [2, n-1). For very small ranges fall back
            # to enumerating the available bases.
            available = n - 3  # count of integers in [2, n-2]
            if available >= k:
                try:
                    bases = random.sample(range(2, n - 1), k)
                except (ValueError, OverflowError):
                    # Fallback for very large `n` where sampling the full range
                    # would overflow internal C ssize_t lengths: generate k
                    # random integers using randrange which supports big ints.
                    bases = [random.randrange(2, n - 1) for _ in range(k)]
            else:
                # small domain: test all possible bases in range
                bases = list(range(2, n - 1))

    # Normalize bases: map large base constants into the valid range [2, n-2]
    # by taking modulus and filter out any values that are not in that range.
    norm_bases = []
    seen = set()
    for b in bases:
        a = b % n
        if a < 2 or a > n - 2:
            continue
        if a in seen:
            continue
        seen.add(a)
        norm_bases.append(a)

    # If normalization removed all bases (possible for very small n),
    # fall back to a small deterministic set within range.
    if not norm_bases:
        fallback = [2, 3, 5, 7, 11]
        norm_bases = [a for a in fallback if a <= n - 2]

    # 4. Perform the test for each base in `norm_bases`.
    for a in norm_bases:
        # Calculate initial term: x = a^d mod n (using fast built-in pow)
        x = power(a, d, n)
        
        # Check Step 2: If x == 1 or x == n - 1, it passes for this base.
        if x == 1 or x == n - 1:
            continue
            
        # Check Step 3: Repeated squaring
        for _ in range(s - 1):
            x = (x * x) % n
            
            # If x becomes n - 1 (which is -1 mod n), it passes for this base.
            if x == n - 1:
                break
                
        # If the inner loop finishes without finding x == n - 1, 
        # then n is definitely composite for this base 'a'.
        if x != n - 1:
            return False # n is a composite number (a is a witness)
            
    # If n passes all k tests, it is highly likely a prime number.
    return True

def is_probable_prime(n, k=8):
    """Convenience wrapper for `is_prime_miller_rabin`.

    Uses deterministic bases automatically for values that fit in 64 bits
    to make the test deterministic on that domain. For larger `n` the test
    remains probabilistic by default.
    """
    if n < (1 << 64):
        return is_prime_miller_rabin(n, k=k, deterministic=True)
    return is_prime_miller_rabin(n, k=k, deterministic=False)

def next_probable_prime(start: int) -> int:
    """
    Finds the next probable prime greater than or equal to 'start'.

    Args:
        start (int): The starting integer to search for the next prime.  

    Returns:
        int: The next probable prime number >= start.
    """
    if start <= 2:
        return 2
    # Ensure we start with an odd number
    candidate = start if start % 2 != 0 else start + 1
    
    while True:
        if is_probable_prime(candidate):
            return candidate
        candidate += 2  # Check only odd numbers


def previous_probable_prime(start: int) -> int:
    """
    Finds the previous probable prime less than or equal to `start`.

    Args:
        start (int): The starting integer to search for the previous prime.

    Returns:
        int: The previous probable prime number <= start. If `start <= 2`
             returns 2.
    """
    if start <= 2:
        return 2

    # Ensure we start with an odd number <= start
    candidate = start if start % 2 != 0 else start - 1

    while candidate >= 3:
        if is_probable_prime(candidate):
            return candidate
        candidate -= 2  # Check only odd numbers going downwards

    # If none found (shouldn't happen), return 2 as a safe fallback
    return 2

## Unit tests
# ----------------------------------------------------------------------------
import unittest


class TestMillerRabin(unittest.TestCase):
    def test_known_primes_and_composites(self):
        p1 = 179424691  # A relatively large prime
        c1 = 179424690  # A composite number
        p2 = 999983     # A smaller prime
        c2 = 561        # Carmichael number (should be composite under Miller-Rabin with sufficient rounds)

        self.assertTrue(is_probable_prime(p1))
        self.assertFalse(is_probable_prime(c1))
        self.assertTrue(is_probable_prime(p2))
        self.assertFalse(is_probable_prime(c2))

        self.assertTrue(is_prime_miller_rabin(7))
        self.assertFalse(is_prime_miller_rabin(9))

    def test_next_and_previous_probable_prime(self):
        # next/previous helpers
        self.assertEqual(next_probable_prime(14), 17)
        self.assertEqual(previous_probable_prime(100), 97)
        self.assertEqual(next_probable_prime(2), 2)
        self.assertEqual(previous_probable_prime(2), 2)

    def test_large_prime_gap(self):
        # Fail immediately unless `high` is the next prime
        # after `low` or `low` is the previous prime before `high`.
        low = 4652353
        high = 4652507

        next_after_low = next_probable_prime(low + 1)
        prev_before_high = previous_probable_prime(high - 1)

        self.assertTrue(
            next_after_low == high or prev_before_high == low,
            msg=(f"Bounds are not adjacent primes: next_after_low={next_after_low}, "
                 f"prev_before_high={prev_before_high}"),
        )


class TestMillerRabinSlow(unittest.TestCase):
    """Slow, expensive tests. Run only when RUN_SLOW=1 is set (or --run-slow passed).
    """

    def test_large_mersenne_127(self):
        # Skip at runtime if slow tests are not enabled. Decorators were
        # previously used but were evaluated at import time (before
        # __main__ could set RUN_SLOW), so do the check inside the test.
        if os.getenv('RUN_SLOW') != '1':
            self.skipTest('skip slow tests by default')
        # 2**127 - 1 is a known Mersenne prime; very expensive to check
        n = 2**127 - 1
        # large k for high confidence; slow by design
        self.assertTrue(is_probable_prime(n, k=200))

    def test_large_k_mersenne(self):
        if os.getenv('RUN_SLOW') != '1':
            self.skipTest('skip slow tests by default')
        # Same candidate with a very large k to stress the algorithm
        n = 2**521 - 1
        self.assertTrue(is_probable_prime(n, k=500))

    def test_carmichael_many_bases(self):
        if os.getenv('RUN_SLOW') != '1':
            self.skipTest('skip slow tests by default')
        # Run Miller-Rabin with many bases on several (composite) Carmichael numbers
        # Expect False for composites even with many bases
        # include several well-known Carmichael numbers (small -> medium-large)
        carmichaels = [561, 1105, 1729, 2465, 2821, 41041, 825265, 321197185, 
                       5394826801, 1436697831295441]
        for c in carmichaels:
            self.assertFalse(is_probable_prime(c, k=500))

# ----------------------------------------------------------------------------
# If run as a script, execute the unit tests
# ----------------------------------------------------------------------------
if __name__ == '__main__':
    import sys
    import argparse
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--run-slow', action='store_true', help='Include slow tests')
    # parse known args and leave the rest for unittest
    args, remaining = parser.parse_known_args()

    # Pass through any remaining args to unittest (verbosity, pattern, etc.)
    # If requested, set RUN_SLOW so slow tests marked with skipUnless run
    if args.run_slow:
        os.environ['RUN_SLOW'] = '1'

    unittest.main(argv=[sys.argv[0]] + remaining)
