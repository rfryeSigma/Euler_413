### THEORY ###
Let us define new variables for equation a^4 + b^4 = d^4 - c^4 .
We can choose c, d odd and a, b = 0 mod 8, with a < b. 
Define:
    t = a / 8, u = b / 8, v = (d - c) / 2, w = (d + c) / 2.
Then the equation becomes
    (t^4 y + u^4) * 2^9 = m = v * w * (v^2 + w^2) = v^3 * w + v * w^3
We search for a candidate (t, u) and calculate the 4th power sum m.
We factor m just enough so that either v or w has a high probability
of being completely factored.
We search for a small subset of the factors with with product x and 
co-product m / x = y. Then we treat x as a constant in the monotonically
increasing cubic function
    f(y) = x * y^3 + x^3 * y - m / x
We set to zero and solve for the root y using binary search or Newton's method.
Knowing x and y, we convert them back to c and d for a solution.

### FOLDERS ###
The code in this folder uses this method to find solutions.

The code in the supplement sub-folder is used to supplement this code. See it's README.

The data files in the several tree folders hold .bin files for product trees.

### CODE ###
Files in this folder: (not completely designed yet)
The main entry point will be search_413.py, which will use a factoring and a solving script.
