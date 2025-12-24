""" 
Select a fast factoring strategy optimized for
    a^4 + b^4 = c^4 + d^4 ~ 10^32 when d +- c ~ 10^8.
to be run on my MAC laptop with M1 chip.
I could use sympy or primefac packages if I wanted full factorization,
but I just want a high probability that one of the algebraic factors
d + c or d - c is completely factored.

I plan to remove common factors of a and b first. These can be any primes.
When factoring the sum of 4th powers, the primes must be odd.

"""
import math
import random
import time
from probable_prime import is_probable_prime