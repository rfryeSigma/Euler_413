from math import gcd, prod
from random import randint

primes_250 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583]
small_limit = 50
small_primes = frozenset(primes_250[:small_limit])
M = prod(small_primes)
large_limit = 150
large_primes = frozenset(primes_250[small_limit:large_limit])
M_large = prod(large_primes)

def miller_rabin(n, bases, pow=pow):
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

def is_prime(n:int, gcd=gcd, randint=randint, 
            miller_rabin=miller_rabin):
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
    if g > 1: return (n == g) and n in small_primes
    g = gcd(n, M_large)
    if g > 1: return (n == g) and n in large_primes

    # 3. Deterministic Range (up to 2^64). First 12 primes ok, but 7 faster
    DETERMINISTIC_BASES_64 = [2, 325, 9375, 28178, 450775, 9780504, 1795265022]
    if n < 18_446_744_073_709_551_616:
        return miller_rabin(n, DETERMINISTIC_BASES_64)
    
    # 4. Extended Range (up to 10^32)
    # We use deterministic bases + 5 random bases for "Computational Certainty"
    extended_bases = DETERMINISTIC_BASES_64 + [randint(2, n-2) for _ in range(5)]
    return miller_rabin(n, extended_bases)


# The positions in a 2-3-5 wheel (period of 30)
# These are all numbers n where gcd(n, 30) == 1
WHEEL_POS = [1, 7, 11, 13, 17, 19, 23, 29]
# Forward gaps between WHEEL_POS plus 31-29
WHEEL_GAPS = [6, 4, 2, 4, 2, 4, 6, 2]
# Reverse gaps from 31 to 29, 29 to 23, ..., 11 to 7, 7 to 1
PREV_GAPS = [2, 6, 4, 2, 4, 2, 4, 6] 
 
def next_prime(start: int) -> int:
    """Finds the next prime >= start using a 2-3-5 wheel."""
    if start > 5:
        # Get onto the wheel
        base = (start // 30) * 30
        rem = start % 30
        
        idx = 0
        for i, p in enumerate(WHEEL_POS):
            if base + p >= start:
                idx = i
                break
        else:
            # If no position in this block works, move to the first position of next block
            base += 30
            idx = 0        
        candidate = base + WHEEL_POS[idx]
        
        # Search forward
        while True:
            if is_prime(candidate): return candidate
            candidate += WHEEL_GAPS[idx]
            idx = (idx + 1) % 8
        assert False

    # Fallback for small numbers
    if start <= 3: return 3 if start == 3 else 2
    return 5

def previous_prime(start: int) -> int:
    """Finds the previous prime <= start using a 2-3-5 wheel."""
    # Small primes handled manually because the wheel skips them
    if start > 7:
        # Valid positions relative to modulo 30:
        # 1, 7, 11, 13, 17, 19, 23, 29
        positions = [1, 7, 11, 13, 17, 19, 23, 29]
        
        # 1. Snap start to the wheel
        base = (start // 30) * 30
        rem = start % 30
        
        # Find the largest index i such that positions[i] <= rem
        idx = -1
        for i in range(len(positions) - 1, -1, -1):
            if positions[i] <= rem:
                idx = i
                break
        
        # If rem < 1, we need to wrap around to the previous block of 30
        if idx == -1:
            base -= 30
            idx = 7 # point to 29
            
        start = base + positions[idx]

        # 2. Iterate backwards
        while start > 7:
            if is_prime(start):
                return start
            
            idx -= 1
            if idx < 0:
                base -= 30
                idx = 7
            start = base + positions[idx]

    # Final fallbacks for primes below the wheel's first entry (7)
    if start >= 5: return 7 if start == 7 else 5
    return 3 if start >= 3 else 2

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
            next_after_low == high and prev_before_high == low,
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
        self.assertTrue(is_prime(n))

    def test_large_k_mersenne(self):
        # Same candidate with a very large k to stress the algorithm
        n = 2**521 - 1
        self.assertTrue(is_prime(n))

    def test_carmichael_many_bases(self):
        # Run Miller-Rabin with many bases on several (composite) Carmichael numbers
        # Expect False for composites even with many bases
        # include several well-known Carmichael numbers (small -> medium-large)
        carmichaels = [561, 1105, 1729, 2465, 2821, 41041, 825265, 321197185, 
                       5394826801, 1436697831295441]
        for c in carmichaels:
            self.assertFalse(is_prime(c))
    
    def test_next_prime(self):
        p = 3
        for _ in range(1_000_000 - 2):
            p = next_prime(p + 1)
        self.assertEqual(p, 15_485_863)
        for _ in range(1_000_000 - 2):
            p = previous_prime(p - 1)
        self.assertEqual(p, 3)

if __name__ == '__main__':
    unittest.main()