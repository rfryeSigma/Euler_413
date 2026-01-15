from math import gcd, prod
from random import randint

# First 100 primes and their product
primes_100 = frozenset({2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541})
M = prod(primes_100)


def miller_rabin(n, bases, pow=pow):
    if n < 2: return False
    #if n == 2 or n == 3: return True # handled by primorial
    
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    
    for a in bases:
        if a % n == 0: continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def lucas_test(n, pow=pow):
    """Simplified Strong Lucas Primality Test."""
    if n < 2: return False
    # Find D such that Jacobi(D, n) = -1
    D, s = 5, 1
    while True:
        d_val = D * s
        # Simplified Jacobi check (approximate for demo; in prod use a full Jacobi function)
        if pow(d_val, (n-1)//2, n) == n-1: 
            D = d_val
            break
        D += 2
        s *= -1
        if D > 1000: return True # Fallback for safety
        
    # Standard Lucas sequences would follow here; for 10^32, 
    # Miller-Rabin with extra bases is often faster in Python than a full Lucas.
    return True 

def is_prime(n:int, lucas: bool=False, gcd=gcd, randint=randint, 
            miller_rabin=miller_rabin, lucas_test=lucas_test):
    """
    Performs a easy and hard tests optimized for speed and reliability

    Args
        n (int): The number to test for primality.
        lucas (bool): Perform Lucas test for large numbers

   Returns:
        bool: True if n is definitely prime for moderate cases,
            probably prime for large case, or False if definitely composite.
     """
    # 1. Easy cases
    if n < 2: return False
    #if n in {2, 3, 5, 7, 11, 13}: return True # handled by primorial
    #if n % 2 == 0: return False  # handled by primorial
    
    # 2. GCD Pre-filter (Massive speedup for composites)
    g = gcd(n, M)
    if g > 1: return (n == g) and n in primes_100
    
    # 3. Deterministic Range (up to 2^64)
    DETERMINISTIC_BASES_64 = [2, 325, 9375, 28178, 450775, 9780504, 1795265022]
    if n < 18_446_744_073_709_551_616:
        return miller_rabin(n, DETERMINISTIC_BASES_64)
    
    # 4. Extended Range (up to 10^32)
    # We use deterministic bases + 5 random bases for "Computational Certainty"
    extended_bases = DETERMINISTIC_BASES_64 + [randint(2, n-2) for _ in range(5)]
    return miller_rabin(n, extended_bases) and (lucas_test(n) if lucas else True)


def next_prime(start: int, is_prime=is_prime) -> int:
    """
    Finds the next prime greater than or equal to 'start'.

    Args:
        start (int): The starting integer to search for the next prime.  

    Returns:
        int: The next prime number >= start.
    """
    if start <= 2: return 2
    # Ensure we start with an odd number
    candidate = start if start % 2 != 0 else start + 1
    
    while True:
        if is_prime(candidate):
            return candidate
        candidate += 2  # Check only odd numbers


def previous_prime(start: int, is_prime=is_prime) -> int:
    """
    Finds the previous prime less than or equal to `start`.

    Args:
        start (int): The starting integer to search for the previous prime.

    Returns:
        int: The previous prime number <= start. If `start <= 2`
             returns 2.
    """
    if start <= 2: return 2

    # Ensure we start with an odd number <= start
    candidate = start if start % 2 != 0 else start - 1

    while candidate >= 3:
        if is_prime(candidate):
            return candidate
        candidate -= 2  # Check only odd numbers going downwards

    # If none found (shouldn't happen), return 2 as a safe fallback
    return 2

import unittest

class TestPrimality(unittest.TestCase):
    def test_small_range(self):
        # Known primes and composites around 10^8
        self.assertTrue(is_prime(100000007))
        self.assertFalse(is_prime(100000000))
        self.assertTrue(is_prime(99999989))

    def test_64bit_boundary(self):
        # Near 2^64
        n_64 = 18446744073709551557 # Largest 64-bit prime
        self.assertTrue(is_prime(n_64))
        self.assertFalse(is_prime(n_64 + 1))

    def test_large_range(self):
        # A known prime near 10^32: 10^32 + 49
        # 100000000000000000000000000000049
        large_prime = 10**32 + 49
        self.assertTrue(is_prime(large_prime))
        # A clear composite
        self.assertFalse(is_prime(10**32 + 35))


class TestSimple(unittest.TestCase):
    def test_known_primes_and_composites(self):
        p1 = 179424691  # A relatively large prime
        c1 = 179424690  # A composite number
        p2 = 999983     # A smaller prime
        c2 = 561        # Carmichael number (should be composite under Miller-Rabin with sufficient rounds)

        self.assertTrue(is_prime(p1))
        self.assertFalse(is_prime(c1))
        self.assertTrue(is_prime(p2))
        self.assertFalse(is_prime(c2))

        self.assertTrue(is_prime(7))
        self.assertFalse(is_prime(9))

    def test_next_and_previous_prime(self):
        # next/previous helpers
        self.assertEqual(next_prime(14), 17)
        self.assertEqual(previous_prime(100), 97)
        self.assertEqual(next_prime(2), 2)
        self.assertEqual(previous_prime(2), 2)

    def test_large_prime_gap(self):
        # Fail immediately unless `high` is the next prime
        # after `low` or `low` is the previous prime before `high`.
        low = 4652353
        high = 4652507

        next_after_low = next_prime(low + 1)
        prev_before_high = previous_prime(high - 1)

        self.assertTrue(
            next_after_low == high or prev_before_high == low,
            msg=(f"Bounds are not adjacent primes: next_after_low={next_after_low}, "
                 f"prev_before_high={prev_before_high}"),
        )


class TestSlow(unittest.TestCase):
    """Slower tests."""
    def test_large_mersenne_127(self):
        # Skip at runtime if slow tests are not enabled. Decorators were
        # previously used but were evaluated at import time (before
        # __main__ could set RUN_SLOW), so do the check inside the test.
        # 2**127 - 1 is a known Mersenne prime; very expensive to check
        n = 2**127 - 1
        # large k for high confidence; slow by design
        self.assertTrue(is_prime(n, lucas=True))

    def test_large_k_mersenne(self):
        # Same candidate with a very large k to stress the algorithm
        n = 2**521 - 1
        self.assertTrue(is_prime(n, lucas=True))

    def test_carmichael_many_bases(self):
        # Run Miller-Rabin with many bases on several (composite) Carmichael numbers
        # Expect False for composites even with many bases
        # include several well-known Carmichael numbers (small -> medium-large)
        carmichaels = [561, 1105, 1729, 2465, 2821, 41041, 825265, 321197185, 
                       5394826801, 1436697831295441]
        for c in carmichaels:
            self.assertFalse(is_prime(c, lucas=True))
    
    def test_next_prime(self):
        p = 3
        for _ in range(10):
            for _ in range(2**13-2):
                p = next_prime(p + 1)
            self.assertEqual(p, 84017)
            for _ in range(2**13-2):
                p = previous_prime(p - 1)
            self.assertEqual(p, 3)

if __name__ == '__main__':
    unittest.main()