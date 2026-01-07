"""
Search for solutions to equation a^4 + b^4 = d^4 - c^4 .
We can choose c, d odd and a, b = 0 mod 8, with a < b. 
Define:
    t = a / 8, u = b / 8, v = (d - c) / 2, w = (d + c) / 2.
Then the equation becomes
    (t^4 y + u^4) * 2^9 = m = v * w * (v^2 + w^2) = v^3 * w + v * w^3
We search for a pair t < u such that one or both is 0 mod 5 but not
both 0 mod 25. The pair must not have any common factors except 
2 and 5 in order to be minimal. We calculate the 4th power sum m.
We factor m just enough so that either v or w has a high probability
of being completely factored.
We search for a small subset of the factors with with product x and 
co-product m / x = y. Then we treat x as a constant in the monotonically
increasing cubic function
    f(y) = x * y^3 + x^3 * y - m / x
We set to zero and solve for the root y using binary search or Newton's method.
Knowing x and y, we convert them back to c and d for a solution.
"""
import concurrent.futures
from multiprocessing import Pool
from time import time
from factoring import factor_tu
from logging_413 import V, IntFlag, parse_flags

def process_tu(t: int, u: int, vV: IntFlag=V.NONE) -> None:
    """Factor and solve for (t, u).
    Will abort with message if find solution"""
    V.log(vV, V.PROCESS, f'Process t {t} in u {u}')
    factoring_results = factor_tu(t, u, vV)
    # CALL SOLVING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    return None

def search_inner_loop(u, vV: IntFlag=V.NONE):
    u_mod_25 = u % 25
    V.log(vV, V.OUTER, f'Outer step u {u} state {u_mod_25}')
        
    # 1. Define the 'Wheel' for t based on u's state
    if u_mod_25 == 0:
        # Rule: Skip t if t % 25 == 0
        t_generator = (t for t in range(1, u) if t % 25 != 0)
        
    elif u_mod_25 in (5, 10, 15, 20):
        # Rule: No restriction on t
        t_generator = range(1, u)
        
    elif u % 5 != 0: # but not 0 mod 25
        # Rule: t % 5 must be 0
        t_generator = range(5, u, 5)
        
    else:
        assert False, "Imposible u case"

    # 2. Execute the search in the selected wheel
    for t in t_generator:
        process_tu(t, u, vV)

def search_p_ut(first_u :int, last_u :int, workers :int=8,
           vV: IntFlag=V.NONE) -> None:
    """
    Performs the search in parallel using a wheel-based lookup for
    t (smallest variable) based on the successive modular states
    of range on u. (second smallest variable)
    """
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        # Map the outer loop values to the executor
        # This returns a dictionary of {future: outer_value}
        future_to_val = {executor.submit(search_inner_loop, u, vV):
                         u for u in range(first_u, last_u + 1)}
        # This loop yields results immediately as they finish
        for future in concurrent.futures.as_completed(future_to_val):
            result = future.result()
    """
    with Pool(processes=workers) as pool:
        async_results = [
            pool.apply_async(search_inner_loop, args=(u, vV))
            for u in range(first_u, last_u + 1)
        ]
        for res in async_results:
            _ = res.get() # This waits for results
    """
    None

def search_ut(first_u :int, last_u :int,
           vV: IntFlag=V.NONE) -> None:
    """
    Performs the search using a wheel-based lookup for
    t (smallest variable) based on the successive modular states
    of range on u. (second smallest variable)
    """
    for u in range(first_u, last_u + 1):
        search_inner_loop(u, vV)

def search_t(u :int, first_t: int=1, last_t :int=-1, 
             vV: IntFlag=V.NONE) -> None:
    """
    Performs the search using a wheel-based lookup for
    t (smallest variable) based on the modular state of 
    given u. (second smallest variable)
    """
    u_mod_25 = u % 25
    if last_t == -1: last_t = u - 1
         
    V.log(vV, V.OUTER, f'Select inner generator for u {u} state {u_mod_25}')
    # 1. Define the 'Wheel' for t based on u's state
    if u_mod_25 == 0:
        # Rule: Skip t if t % 25 == 0
        t_generator = (t for t in range(first_t, last_t+1) if t % 25 != 0)
        
    elif u_mod_25 in (5, 10, 15, 20):
        # Rule: No restriction on t
        t_generator = range(first_t, last_t+1)
        
    elif u % 5 != 0: # but not 0 mod 25
        # Rule: t % 5 must be 0
        t_mod_5 = first_t % 5
        first_t += 0 if t_mod_5 == 0 else 5 - t_mod_5
        t_generator = range(first_t, last_t+1, 5)
        
    else:
        assert False, "Imposible u case"

    # 2. Execute the search in the selected wheel
    for t in t_generator:
        process_tu(t, u, vV)

## Unit tests
# ----------------------------------------------------------------------------
import unittest

class TestSearch412Fast(unittest.TestCase):

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

    # search_ut_p command
    p_search_p_ut = subparsers.add_parser('search_p_ut', 
                        help='Search in parallel all t for range on u')
    p_search_p_ut.add_argument('first_u', type=str, 
                        help='Python expression for first in outer range')
    p_search_p_ut.add_argument('last_u', type=str, 
                        help='Python expression for last in outer range')
    p_search_p_ut.add_argument('workers', nargs='?', type=int, default=8,
                        help='Number of parallel workers')
    p_search_p_ut.add_argument('-v', '--vV', type=parse_flags, default=V.NONE,
        help="Set verbosity levels (e.g., 'OUTER', 'OUTER|PROCESS', 'DEBUG')")

    # search_ut command
    p_search_ut = subparsers.add_parser('search_ut', help='Search all t for range on u')
    p_search_ut.add_argument('first_u', type=str, 
                        help='Python expression for first in outer range')
    p_search_ut.add_argument('last_u', type=str, 
                        help='Python expression for last in outer range')
    p_search_ut.add_argument('-v', '--vV', type=parse_flags, default=V.NONE,
        help="Set verbosity levels (e.g., 'OUTER', 'OUTER|PROCESS', 'DEBUG')")

    # search_t command
    p_search_t = subparsers.add_parser('search_t', help='Search range on t for one u')
    p_search_t.add_argument('u', type=str, help='Python expression for single u')
    p_search_t.add_argument('first_t', type=str,
                        help='Python expression for first in inner range')
    p_search_t.add_argument('last_t', type=str, nargs='?', default='-1',
                        help='Python expression for last in inner range')
    p_search_t.add_argument('-v', '--vV', type=parse_flags, default=V.NONE,
        help="Set vV levels (e.g., 'OUTER', 'OUTER|PROCESS', 'DEBUG')")


    """ parse known args and leave the rest for unittest
        -v	--verbose	Show individual test names and results.
        -q	--quiet	Show minimal output.
        -f	--failfast	Stop on first failure.
        -b	--buffer	Hide stdout/stderr for passing tests.
        -c	--catch	Graceful interrupt with Ctrl-C.
        -k	(None)	Run tests matching a pattern
    """
    args, rest = parser.parse_known_args(argv)

    if args.command == 'search_p_ut':
        first_u = eval(args.first_u)
        last_u = eval(args.last_u)
        start = time()
        search_p_ut(first_u, last_u, workers=args.workers, vV=args.vV)
        elapsed = time() - start
        print(f'Elapsed: {elapsed:.6f}s')
        return

    if args.command == 'search_ut':
        first_u = eval(args.first_u)
        last_u = eval(args.last_u)
        start = time()
        search_ut(first_u, last_u, vV=args.vV)
        elapsed = time() - start
        print(f'Elapsed: {elapsed:.6f}s')
        return

    if args.command == 'search_t':
        u = eval(args.u)
        first_t = eval(args.first_t)
        last_t = eval(args.last_t)
        start = time()
        search_t(u, first_t, last_t, vV=args.vV)
        elapsed = time() - start
        print(f'Elapsed: {elapsed:.6f}s')
        return

    # No command, so run unit tests
    return unittest.main(argv=[sys.argv[0]] + rest)

if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
