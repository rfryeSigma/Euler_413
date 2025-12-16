# Euler_413
Solve a^4 + b^4 + c^4 = d^4

This repository contains information on Euler's Sum of Powers problem
for the (4,1,3) pattern. This pattern is the equation with exponent 4,
1 variable on one side, and 3 variables on the other side.

In 1769, Leonhard Euler extended Fermat's Last Theorem by claiming that
for exponent p>2, the sum of p-1 integers raised to the p-th power cannot
equal an integer raised to the p-th power. In 1966, L. J. Lander and
T. R. Parkin found a counter-example for the 5-th power.

In 1988, Noam Elkies proved that there are an infinite number of 
counter-examples by inventing a parametric method combined with Elliptic
Curve search. He illustrated two solutions. In the same year, I found the
minimum counter-example as ranked first by exponent and then by sum.

Since then, 7 other solutions have been catalogued at euler.free.fr .
That web site also catalogues solutions to other sum patterns defined
by a more general 1966 conjecture of Lander, Parkin, and Selfridge.

Aside from Elkies' parametric approach, all brute force searches for
counter-examples start by transposing the c term to the other side,
searching for (d, c) pairs, and exploiting the rich algebraic and modular
restrictions when the difference has to equal a sum on the (a, b) pair.

In my code here, I implement the opposite approach. I search for (a, b)
pairs and thereby give up the hope for many restrictions on candidates.
Then I partially factor the sums and search for a partition of the small
factors to be used as a constant in a monotonically increasing cubic
equation. I solve for the root of the equation and use it to construct
the solution.
