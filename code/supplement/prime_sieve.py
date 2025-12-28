# Hand code fast single buffer segmented sieve of Eratosthenes.

import math
import numpy as np
import os
import time

class SegmentedSieve:
    def __init__(self):
        # Hidden internal buffers
        self._prime_buffer = None
        
        # Publicly accessible
        self.primes = np.empty(0) # view of prime buffer containing primes
        self.extra = np.empty(0) # view of prime buffer after primes
        self.n_primes = 0
        self.primes_1mod8 = np.empty(0) # extract of primes: p%8 == 1

    def run_sieve(self, n_primes: int=1_000_000, overhead: int=-1, 
                  verbose: bool=False) -> int|None:
        """
        Initializes buffers and fills prime_buffer by sieving segments.

        Args:
            n_primes (int): Number of primes to generate.
            overhead (int): Number of extra numbers to allocate in buffer.
                Use -1 to default to min(n_primes, 2**16).

        Returns int | None:
            Value of last prime generated
        """
        start_time = time.time()
        if n_primes == 0:
            return None
        n = self.n_primes = int(n_primes)
        overhead = int(overhead)
        if overhead == -1: overhead = min(n_primes, 65_536)
        capacity = n + int(overhead)
        if verbose:
            print(f"Generating {n:_} primes in {capacity:_} length buffer...")

        # Initialize internal state. unint could handle large values, but int faster.
        pb = self._prime_buffer = np.empty(capacity, dtype=np.int64)
        self.primes = self._prime_buffer[:n]
        self.extra = self._prime_buffer[n:]

        pb[0] = 2

        # Sieve for initial sieving primes in rest of buffer
        is_prime = pb[1:]
        is_prime[:] = 1 # Initially all odd values are candidates.
        # The list of odd candidates is offset so first candidate is the prime 3.
        # So the candidate at index 0 is: k = 2 * i + 3, which can be inverted for
        # the index of a candidate: i = (k - 3) // 2.
        # When filtering multiples of a given prime, the first candidate is p * p
        # with location i = (p * p - 3) // 2. We want the final location sieved to be
        # i <= len(is_prime) - 1. So we can solve the relation:
        # (p_max * p_max -3) // 2 <= len(is_prime) - 1
        p_max = math.isqrt(2 * (len(is_prime) - 1) + 3)
        i_max = (p_max - 3) // 2
        for i in range(i_max + 1):
            if is_prime[i]:
                p = 2 * i + 3
                is_prime[(p * p - 3) // 2 :: p] = 0
        indices = np.nonzero(is_prime)[0]
        if verbose:
            elapsed = time.time() - start_time
            print(f'Pre-sieve found {len(indices):_} odd primes filtering thru {p_max}'
                  f' elapsed={elapsed:.5f}s')
        c = 1 + len(indices)
        pb[1:c] = 2 * indices + 3

        seg_count = 0
        while c < n:
            # 1. Setup offset and segment
            is_prime = pb[c:] # View of the remainder of the buffer
            len_is = len(is_prime)
            is_prime[:] = 1   # Reset segment to 'True'
            
            # 2. Sieve with ALREADY found primes (Vectorized start calculation)
            # We slice up to c to use primes we've already solidified
            sieving_primes = pb[1:c]
            last_sieving_prime = int(pb[c-1])
            offset = last_sieving_prime + (2 if c > 1 else 1)
            p_max = math.isqrt(2 * (len(is_prime) - 1) + offset)
            idx_limit = np.searchsorted(sieving_primes, p_max, side='right')
            
            for p in sieving_primes[:idx_limit]:
                p_val = int(p)
                # Faster modular math for the first index
                p0 = (-offset % p_val)
                if p0 % 2 != 0: p0 += p_val
                p0 //= 2
                if p0 < len_is:
                    is_prime[p0::p_val] = 0

            # 3. Finalize the new primes
            # The nonzero indices now represent confirmed primes
            new_prime_indices = np.nonzero(is_prime)[0]
            new_primes = 2 * new_prime_indices + offset
            
            if len(new_primes) == 0:
                print(f'ERROR only {c}/{n} primes with overhead {overhead:_}')
                break
                
            num_new = min(len(new_primes), n - c)
            pb[c : c + num_new] = new_primes[:num_new]
            c += num_new
            seg_count += 1
            
            if verbose:
                elapsed = time.time() - start_time
                print(f"[sieve] segment {seg_count} "
                      f"primes={c:_}/{n:_} "
                      f"({(100 * c / n):.1f}%) elapsed={elapsed:.5f}s")

        return self.primes[-1]

    def select_primes_1mod8(self) -> int:
        """
        Select primes from prime_buffer such that p%8 == 1. Save in object.

        Returns int
            Number of primes selected
        """
        self.primes_1mod8 = [p for p in self.primes if p%8 == 1]
        return len(self.primes_1mod8)

### END class SegmentedSieve


class PrimeProductTree:
    @staticmethod
    def build(primes: list)-> list:
        """Constructs a binary product tree from a list of primes.
        The tree has the primes as the leaves with each layer
        a list of the pairwise product of previus layer.
        The root is the product of all.

        Args:
            primes (list): List of primes for use as leaves of tree.
                Most efficient if length is a power of 2.
                Else passes through odd remainders.

        Returns list:
            List of lists of products
        """
        len_p = len(primes)
        assert 0 < len_p
        if len_p.bit_count() != 1:
            print(f'WARNING: number of primes {len_p} is not a power of 2')
        len_t = len_p.bit_length() + (len_p.bit_count() > 1)
        tree = [[]] * len_t
        lvl = 0
        tree[lvl] = curr = primes
        while len(curr) > 1:
            next_lvl= [0] * math.ceil(len(curr) / 2)
            for i in range(len(next_lvl)):
                i2 = 2 * i
                next_lvl[i] = curr[i2] * (curr[i2+1] if i2 + 1 < len(curr) else 1)
            lvl += 1
            tree[lvl] = curr = next_lvl
        return tree

    @staticmethod
    def save(tree: list, folder: str) -> None:
        """
        Save a product tree of primes as a folder containing a file for each layer.
        
        Args:
            tree (list): prime product tree to write to disk.
            folder (str): name of folder to hold the prime tree layer files.
        """
        if not os.path.exists(folder):
            os.makedirs(folder)

        # Save leaf primes (level 0) and internal levels
        for i in range(len(tree)):
            with open(f"{folder}/level_{i}.bin", "wb") as f:
                for node in tree[i]:
                    node_bytes = node.to_bytes((node.bit_length() + 7) // 8, 'little')
                    f.write(len(node_bytes).to_bytes(4, 'little'))
                    f.write(node_bytes)
        
        # Save a metadata file so we know how many levels there are
        with open(f"{folder}/meta.txt", "w") as f:
            f.write(str(len(tree)))


    @staticmethod
    def load(level: int, folder: str) -> list:
        """
        Load all or a single level of a prime product tree.
        The 0th level is the leaves, the primes themselves.
        Intermediate levels are pairwise products of previous level
        The last level is a list containing the product of all primes.
        
        Args:
            level (int): level to load;
                if -1, load root;
                if -2 load all levels
            folder (str): name of folder containing level files
        
        Returns: 
            list of products at selected level or list of lists.
        """
        def load_level(level_idx: int) -> list:
            """Loads a specific level into memory."""
            num_nodes = 2** (num_levels - level_idx - 1)
            nodes = [[]] * num_nodes
            path = f"{folder}/level_{level_idx}.bin"
            with open(path, "rb") as f:
                for node_idx in range(num_nodes):
                    len_bytes = f.read(4)
                    if not len_bytes:
                        return nodes[:node_idx]
                        break
                    length = int.from_bytes(len_bytes, 'little')
                    nodes[node_idx] = (int.from_bytes(f.read(length), 'little'))
                else:
                    assert node_idx + 1 == num_nodes
                return nodes

        def get_root():
            """Returns the single product at the top of the tree."""
            return load_level(num_levels - 1)[0]

        with open(f"{folder}/meta.txt", "r") as f:
            num_levels = int(f.read())
            print(f'{folder} has {num_levels} levels')

            if 0 <= level:
                nodes = load_level(level)
                return nodes
            
            if -1 == level:
                return get_root()
            
            all_nodes = [[]] * num_levels
            for level in range(num_levels):
                all_nodes[level] = load_level(level)
            return all_nodes

    @classmethod
    def factor(cls, target_gcd: int, tree: list, level: int|None=None, index: int=0) -> list:
        """
        Factor target_gcd by recursively descending prime product tree.
        """
        if level is None:
            level = len(tree) - 1 # Start at the top (the root product)

        # Base case: we reached the leaves of the tree (the primes themselves)
        if level == 0:
            p = tree[0][index]
            return [p] if target_gcd % p == 0 else []

        found_primes = []
        
        # Check the "Left" child
        left_index = index * 2
        left_prod = tree[level-1][left_index]
        # Only descend if the GCD is greater than 1
        if math.gcd(target_gcd, left_prod) > 1:
            found_primes.extend(cls.factor(target_gcd, tree, level-1, left_index))

        # Check the "Right" child
        right_index = index * 2 + 1
        if right_index < len(tree[level-1]): # check pass through if level not power of 2
            right_prod = tree[level-1][right_index]
            if math.gcd(target_gcd, right_prod) > 1:
                found_primes.extend(cls.factor(target_gcd, tree, level-1, right_index))

        return found_primes

### END class PrimeProductTree

## Unit tests
# ----------------------------------------------------------------------------
import unittest

class TestSieveFast(unittest.TestCase):

    def test_overhead_fast(self):
        """Check how much overhead needed for n primes
        """
        def check(sieve, n, over, prime):
            p = sieve.run_sieve(n, over)
            self.assertEqual(p, prime)
            self.assertEqual(n, len(sieve.primes))
            
        sieve = SegmentedSieve()

        over = 0 # overhead buffer length
        check(sieve, 1, over, 2) # initial prime
        check(sieve, 2, over, 3) # no prime gap between 2, 3
        check(sieve, 3, over, 5) # no odd prime gap between 3, 5
        check(sieve, 4, over, 7) # no odd prime gap between 5, 7
        p = sieve.run_sieve(5, over) # overhead too small
        self.assertNotEqual(p, 11)

        over = 1
        check(sieve, 4, over, 7) # 4th prime still 7 with over 1
        check(sieve, 5, over, 11) # odd prime gap between 7, 11
        check(sieve, 6, over, 13) # no extra odd prime gap between 11, 13
        check(sieve, 7, over, 17) # odd prime gap between 13, 17 
        check(sieve, 8, over, 19) # no extra odd prime gap between 17, 19
        check(sieve, 9, over, 23) # odd prime gap between 19, 23 
        p = sieve.run_sieve(10, over) # overhead too small
        self.assertNotEqual(p, 29)

        over = 2
        check(sieve, 4, over, 7) # 4th prime still 7 with over 2
        check(sieve, 10, over, 29) # odd prime gap 3 between 23, 29
        check(sieve, 11, over, 31)
        check(sieve, 12, over, 37)
        check(sieve, 13, over, 41)
        check(sieve, 14, over, 43)
        check(sieve, 15, over, 47)
        check(sieve, 16, over, 53)
        check(sieve, 17, over, 59)
        check(sieve, 18, over, 61)
        check(sieve, 19, over, 67)
        check(sieve, 20, over, 71)
        check(sieve, 21, over, 73)
        check(sieve, 22, over, 79)
        check(sieve, 23, over, 83)
        check(sieve, 24, over, 89)
        p = sieve.run_sieve(25, over) # overhead too small
        self.assertNotEqual(p, 97)

        over = 3
        check(sieve, 25, over, 97)

    def test_25_fast(self):
        """Reheck first 25 primes with ample overhead
        """
        def check(sieve, n, over, prime):
            p = sieve.run_sieve(n, over)
            self.assertEqual(p, prime)
            self.assertEqual(n, len(sieve.primes))
            
        sieve = SegmentedSieve()
        over = 3 # overhead needed for 25th prime 97
        primes = (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97)
        self.assertEqual(len(primes), 25)
        for i, p in enumerate(primes, start=1):
            check(sieve, i, over, p)

    def test_twos_fast(self):
        """Test 2**5, 2**10, 2**15, 2**20 primes"""
        pairs = ((2**5, 131), (2**10, 8161), (2**15, 386_093),
                (2**20, 16_290_047))
        sieve = SegmentedSieve()
        for (n, p) in pairs:
            prime = sieve.run_sieve(n)
            self.assertEqual(p, prime)

    def test_tens_fast(self):
        """Test 100, 1000, ..., 1_000_000 primes"""
        pairs = ((100, 541), (1000, 7919), (10_000, 104_729), 
                 (100_000, 1_299_709), (1_000_000, 15_485_863),
                 )
        sieve = SegmentedSieve()
        for (n, p) in pairs:
            prime = sieve.run_sieve(n)
            self.assertEqual(p, prime)

class TestSieveSlow(unittest.TestCase):

    def test_twos_slow(self):
        """Test 2**25, 2**27 primes"""
        pairs = ((2**25, 645_155_197), # takes 2.3s
                #(2**27, 2_777_105_129), # takes 16.3s
                )
        sieve = SegmentedSieve()
        for (n, p) in pairs:
            prime = sieve.run_sieve(n)
            self.assertEqual(p, prime)

    def test_tens_slow(self):
        """Test 10_000_000, ... primes"""
        pairs = ((10_000_000, 179_424_673), # takes 0.5s
                 #(100_000_000, 2_038_074_743), # takes 9.9s
                 )
        sieve = SegmentedSieve()
        for (n, p) in pairs:
            prime = sieve.run_sieve(n)
            self.assertEqual(p, prime)

class Test1Mod8_fast(unittest.TestCase):

    def test_1mod8_fast(self):
        """Test select_primes_1mod8"""
        sieve = SegmentedSieve()
        _ = sieve.run_sieve(12_345)
        c = sieve.select_primes_1mod8()
        self.assertEqual(c, 3050)
 
class Test1Tree_fast(unittest.TestCase):

    def test_build_save_load_fast(self):
        """Test build, save and load prime product tree"""
        tree = PrimeProductTree.build([2,3,5,7,11,13,17])
        self.assertEqual(tree[-1], [510510])
        self.assertEqual(tree[-2], [210, 2431])
        PrimeProductTree.save(tree, 'test')
        t2 = PrimeProductTree.load(-2, 'test')
        self.assertEqual(t2, tree)
        g = PrimeProductTree.load(-1, 'test')
        self.assertEqual(g, 510510)

    def test_factor_fast(self):
        """Test factoring with gcd in prime product tree"""
        tree = tree = PrimeProductTree.build([2,3,5,7,11,13,17])
        factors = PrimeProductTree.factor(19, tree)
        self.assertEqual(factors, [])
        factors = PrimeProductTree.factor(77, tree)
        self.assertEqual(factors, [7, 11])
        factors = PrimeProductTree.factor(17*17, tree)
        self.assertEqual(factors, [17])
        factors = PrimeProductTree.factor(3*17, tree)
        self.assertEqual(factors, [3, 17])

### END unittests 

## --- Main Section ---
import sys
import argparse
def main(argv=None):
    """Command-line dispatcher"""
    parser = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers(dest='command')

    # sieve command
    p_sieve = subparsers.add_parser('sieve', help='Sieve for primes')
    p_sieve.add_argument('n_primes', type=str, 
                        help='Python expression for integer nmber of primes')
    p_sieve.add_argument('overhead', nargs='?', type=str, default='-1',
                        help='Python expression for integer overhead')
    p_sieve.add_argument('-v', action='store_true', help='whether to log progress')

    # p_1mod8 command
    p_p_1mod8 = subparsers.add_parser('p_1mod8', help='Select primes: p%8 == 1')

    """ parse known args and leave the rest for unittest
        -v	--verbose	Show individual test names and results.
        -q	--quiet	Show minimal output.
        -f	--failfast	Stop on first failure.
        -b	--buffer	Hide stdout/stderr for passing tests.
        -c	--catch	Graceful interrupt with Ctrl-C.
        -k	(None)	Run tests matching a pattern
    """
    args, rest = parser.parse_known_args(argv)
    if args.command == 'sieve':
        n_primes = eval(args.n_primes)
        overhead = eval(args.overhead)
        sieve = SegmentedSieve()
        prime = sieve.run_sieve(n_primes, overhead, verbose=args.v)
        print(f'last prime {prime:_}')
        if n_primes < 15:
            print('\t', sieve.primes)
        sys.exit(0)

    if args.command == 'p_1mod8':
        sieve = SegmentedSieve()
        _ = sieve.run_sieve(12_345)
        c = sieve.select_primes_1mod8()
        print(f'Selected {c} 1mod8 primes')
        sys.exit(0)

    # No command, so run unit tests
    return unittest.main(argv=[sys.argv[0]] + rest)

if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
