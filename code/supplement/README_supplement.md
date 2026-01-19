The code and daata in this supplement folder are used to aid the main code.


decompose_known_solutions.py
    Documents the 11 known solutions to the Euler 413 problem.
    Contains code to factor the solutions and look for patterns.

factoring_strategy.py
    Previous version of the main factoring script.
    
prime_sieve.py
    Generates lists of primes and product trees of those primes.

search_413.log
    Hand annotated log of search ranges and results confirmed or found.

PROCESS FOR BUILDING PRIME FACTOR TREES
code % python
>>> import supplement.prime_sieve
>>> sieve = supplement.prime_sieve.SegmentedSieve()
>>> sieve.run_sieve(1e6) # Generate 1_000_000 primes
np.int64(15485863). # returns largest prime in list
        # note the type is not int; doesn't support % operator
>>> sieve.select_primes_mod8()
(249796, 750203). # returns numbers of 1mod8 and 357 primes; saves in sieve
>>> sieve.primes_1mod8[2**16 - 1]
3684641 # last prime in a 17 level (0 thrru 16) prime factor tree
>>> sieve.primes_1mod8[2**17 - 1]
7766713 # last prime in an 18 level (0 thrru 17) prime factor tree
>>> primes = [p for p in map(int, sieve.primes_1mod8[:2**17])]
        # note the map from np.int64 to int is necessary to use % operatore 
>>> len(primes)
131072
>>> primes[0]
17
>>> primes[-1]
7766713
>>> tree = supplement.prime_sieve.PrimeProductTree.build(primes)
>>> len(tree)
18
>>> supplement.prime_sieve.PrimeProductTree.save(tree, 'tree1_18')
>>> 
