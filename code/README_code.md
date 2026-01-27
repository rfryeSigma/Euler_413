### THEORY ###
The goal is to find integer solutions to the equation
    a^4 + b^4 + c^4 = d^4 .

Move the c term to the right side and factor
    d^4 - c^4 = (d - c) * (d + c) * (d^2 + c^2) .

Modulo 8 and accounting for symmetries, the original equation
has only 1 non-trivial solution: 0 + 0 + 1 = 1 ,
so we can choose c, d odd and a, b = 0 mod 8, with a < b. 

Define new variables:
    t = a / 8, u = b / 8, v = (d - c) / 2, w = (d + c) / 2
so that
    a = 2^3 * t and b = 2^3 * u .
And since 2 * v + 2 * w = 2 * d, and 2 * v - 2 * w = 2 * c ,
    d = w + v and c = w - v

Then the left side of the equation becomes
    (2^3 * t)^4 + (2^3 * u)^4 = 2^12 * (t^4 + u^4) = 2^12 * m
And the factored right side becomes
    (2 * v) * (2 * w) * ( (v + w)^2 + (v - w)^2 )
    = 2^3 * v * w * (v^2 + w^2)
    = 2^3 * (v^3 * w + v * w^3)
Notice the symmetry in the final term. If we knew either v or w,
we could call it x and the other variable y.
Then the final term becomes a cubic in y with parameter x:
    x + y^3 + x^2 * y

Putting this all together, we want to solve the equations
    t^4 + u^4 = m; 2^9 * m = x * (y^3 + x^2 * y)

### PROCESS ###
We begin our solution process by searching for pairs t < u .

Modulo 5 and accounting for symmetries, the original equation
has only 1 non-trivial solution: 0 + 0 + 1 = 1 . This is the
same as the modulo 8 solution which allowed us to reduce (a,b) to (t,u).
But the modulo 5 condition is independent, so one or both of (t, u)
must be 0 mod 5. If both t and u were 0 mod 25, the other terms
would have to be the same mod 25, and the solution would be trivial.
These conditions restrict the choice of (t, u).

We factor m just enough so that either v or w has a high probability
of being completely factored.

We select small subsets of the identified factors with product x .
If we can solve the monotonically increasing cubic function
    f(y; x) = y^3 + x^2 * y - 2^9 * m / x
for the root y, then we can convert x and y back to c and d.

The flow of control is that search_413.py selects (t, u) pairs
and calls factoring and then solving to  find (v, w).
When a solution is found, search_414.py raises an assertion error
to interrupt the search process.

### FOLDERS ###
The code in this folder uses this method to find solutions by factoring.

The code in the elliptical sub_folder attempts to match factoring solutions 
with Elliptical Curve solutions. See it's README.

There is a play folder, but it is not registered in GitHub. 
It contains experimental code.

The code in the supplement sub-folder is used to supplement this code. 
See it's README.

The data files in the several tree folders hold .bin files for prime 
product trees.
I am reluctant to include in git because of size,
See instructions for generating in supplement/README_supplement.md

### CODE ###
Files in this folder:

factoring.py:
    Attempts to factor t^4 + u^4 into primes.

logging_413.py:
    Common logging class for multi-case logging. Supports flag parsing.

search_413.py:
    Main entry point. Searches ranges of (t,u). Calls factor and solve.

solutions.py:
    Database of known solutions.

solving.py:
    Partitions factors and attempts to solve cubic equation.

utilities_413.py:
    Tree utilities extracted from supplement/prinme_sieve.py for factoring.
