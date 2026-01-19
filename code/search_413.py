"""
As explained in README_code.md, we search for small pairs t < u such that
one or both is 0 mod 5 but not both 0 mod 25.
For each (t, u) pair, we call the factoring module to partially factor
    m = t^4 + u^4,
and we call the solving module to search for solutions (x, y) to
    m = 2^-9 * (x + y^3 + x^2 * y)
"""
import concurrent.futures
from multiprocessing import Pool
from time import time
from factoring import factor_tu, return_saved_times
from logging_413 import V, IntFlag, parse_flags
from os import system
from solving import solve_factors

def process_tu(t: int, u: int, vV: IntFlag=V.NONE,
               log=V.log, factor_tu=factor_tu, solve_factors=solve_factors, system=system,
               ) -> None:
    """Factor and solve for (t, u).
    Will abort with message if find solution"""
    log(vV, V.PROCESS, f'Process t {t:_} in u {u:_}')
    factoring_results = factor_tu(t, u, vV)
    solving_results = solve_factors(*factoring_results, vV=vV)
    if solving_results is None: return
    system('say found solution')
    assert solving_results is None, \
        f'Found solution (t, u, v, w): ({t:_}, {u:_},' \
        f'{solving_results[0]}, {solving_results[1]})'

def search_inner_loop(u, vV: IntFlag=V.NONE,
                      log=V.log, process_tu=process_tu):
    u_mod_25 = u % 25
    log(vV, V.OUTER, f'Outer step u {u:_} state {u_mod_25}')
        
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

def search_stride_worker(worker_id: int, first_u: int, last_u: int, 
                         workers: int, vV: IntFlag):
    """Processes every Nth 'u'"""
    # The 'stride' is the number of workers
    for u in range(first_u + worker_id, last_u + 1, workers):
        search_inner_loop(u, vV)
    vV.log(vV, V.F_DIAG, f'Worker {worker_id}: {return_saved_times()}')

def search_p_ut(first_u :int, last_u :int, workers :int=8, vV: IntFlag=V.NONE,
           search_stride_worker=search_stride_worker) -> None:
    """
    Performs the search in parallel using a wheel-based lookup for
    t (smallest variable) based on the successive modular states
    of range on u. (second smallest variable)
    """
    with Pool(processes=workers) as pool:
        async_results = [
            pool.apply_async(search_stride_worker, 
                             args=(i, first_u, last_u, workers, vV))
            for i in range(workers)
        ]
        for res in async_results:
            _ = res.get() # This waits for results

def search_ut(first_u :int, last_u :int, vV: IntFlag=V.NONE, 
              search_inner_loop=search_inner_loop) -> None:
    """
    Performs the search using a wheel-based lookup for
    t (smallest variable) based on the successive modular states
    of range on u. (second smallest variable)
    """
    for u in range(first_u, last_u + 1):
        search_inner_loop(u, vV)
    vV.log(vV, V.F_DIAG, return_saved_times())

def search_t(u :int, first_t: int=1, last_t :int=-1, vV: IntFlag=V.NONE,
             log=V.log, process_tu=process_tu) -> None:
    """
    Performs the search using a wheel-based lookup for
    t (smallest variable) based on the modular state of 
    given u. (second smallest variable)
    """
    u_mod_25 = u % 25
    if last_t == -1: last_t = u - 1
         
    log(vV, V.OUTER, f'Select inner generator for u {u:_} state {u_mod_25}')
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
    log(vV, V.F_DIAG, return_saved_times())

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
        system('say Calculation finished')
        return

    if args.command == 'search_ut':
        first_u = eval(args.first_u)
        last_u = eval(args.last_u)
        start = time()
        search_ut(first_u, last_u, vV=args.vV)
        elapsed = time() - start
        print(f'Elapsed: {elapsed:.6f}s')
        system('say Calculation finished')
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
