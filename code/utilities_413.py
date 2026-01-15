# Utility functions used during the 413 process

from math import gcd

def load_tree(level: int, folder: str) -> list:
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
        #print(f'{folder} has {num_levels} levels')

        if 0 <= level:
            nodes = load_level(level)
            return nodes
        
        if -1 == level:
            return get_root()
        
        all_nodes = [[]] * num_levels
        for level in range(num_levels):
            all_nodes[level] = load_level(level)
        return all_nodes

def factor_tree_gen(g: int, tree: list, level: int | None = None, index: int = 0,
                    _func=None, len=len, gcd=gcd):
    """
    Yield prime factors of g by recursively descending prime product tree.
    """
    # Base case: we reached the leaves
    if level == 0:
        p = tree[0][index]
        if g % p == 0: yield p
        return

    # First time
    if level is None:
        level = len(tree) - 1
        _func = factor_tree_gen

    # Check the "Left" child
    left_index = index * 2
    left_prod = tree[level-1][left_index]
    g_left = gcd(g, left_prod)
    if g_left > 1:
        yield from _func(g_left, tree, level-1, left_index, _func=_func)

    # Check the "Right" child
    right_index = index * 2 + 1
    if right_index < len(tree[level-1]):
        right_prod = tree[level-1][right_index]
        g_right = gcd(g, right_prod)
        if g_right > 1:
            yield from _func(g_right, tree, level-1, right_index, _func=_func)

def factor_tree(g: int, tree: list, level: int|None=None, index: int=0,
                _func=None, len=len, gcd=gcd) -> list:
    """
    Factor target_gcd by recursively descending prime product tree.
    """ 
    # Base case: we reached the leaves of the tree (the primes themselves)
    if level == 0:
        p = tree[0][index]
        return [p] if g % p == 0 else []

    # First time
    if level is None:
        level = len(tree) - 1 # Start at the top (the root product)
        _func = factor_tree

    found_primes = []
    
    # Check the "Left" child
    left_index = index * 2
    left_prod = tree[level-1][left_index]
    g_left = gcd(g, left_prod)
    if g_left > 1:
        found_primes.extend(_func(g_left, tree, level-1, left_index, 
                                  _func=_func))

    # Check the "Right" child
    right_index = index * 2 + 1
    if right_index < len(tree[level-1]): # check pass through
        right_prod = tree[level-1][right_index]
        g_right = gcd(g, right_prod)
        if g_right > 1:
            found_primes.extend(_func(g_right, tree, level-1, right_index, 
                                      _func=_func))

    return found_primes

## Unit tests
# ----------------------------------------------------------------------------
import unittest
from math import prod

class TestUtilitiesFast(unittest.TestCase):

    def test_factor_tree_fast(self):
        tree357_13 = load_tree(-2, 'tree357_13')
        primes357 = [421, 115_223]
        prime1 = 17
        n = prod(primes357) * prime1
        result = factor_tree(n, tree357_13)
        self.assertEqual(result, primes357)
        cofactor = n / prod(result)
        self.assertEqual(cofactor, float(prime1))

### END unittests 

## --- Main Section ---
import argparse
import sys
import timeit

def main(argv=None):
    """Command-line dispatcher"""
    parser = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers(dest='command')

    # factor command
    p_factor = subparsers.add_parser('factor', help='Factor with factor_tree')
    p_factor.add_argument('n', type=str, nargs='?',
        default='17*41*73*89*97*193*233*401*409*857*881*1801*1873*3889*929*8377' \
        '*8513*18_233*18_257*39_769*39_841*84_913*84_961*181_777*181_873*389_057' \
        '*389_089*823_777*823_841*1_747_489*1_747_513*3_684_641',
        help='Python expression for number to factor')

    """ parse known args and leave the rest for unittest
        -v	--verbose	Show individual test names and results.
        -q	--quiet	Show minimal output.
        -f	--failfast	Stop on first failure.
        -b	--buffer	Hide stdout/stderr for passing tests.
        -c	--catch	Graceful interrupt with Ctrl-C.
        -k	(None)	Run tests matching a pattern
    """
    args, rest = parser.parse_known_args(argv)

    if args.command == 'factor':
        n = eval(args.n)
        tree1_17 = load_tree(-2, 'tree1_17') # 1mod8, 17 levels, length 2**16, 17 thru 3_684_641

        results = factor_tree(n, tree1_17)
        elapsed = timeit.timeit(lambda: factor_tree(n, tree1_17), 
                                number=1000) / 1000
        cofactor = n / prod(results)
        print('factor_tree:')
        print(f'found {results}\nCofactor {cofactor}\nElapsed: {elapsed:.8f}s')

        results = list(factor_tree_gen(n, tree1_17))
        elapsed = timeit.timeit(lambda: list(factor_tree_gen(n, tree1_17)), 
                                number=1000) / 1000
        cofactor = n / prod(results)
        print('/nfactor_tree_gen:')
        print(f'found {results}\nCofactor {cofactor}\nElapsed: {elapsed:.8f}s')
        return

    # No command, so run unit tests
    return unittest.main(argv=[sys.argv[0]] + rest)

if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
