These are notes about how Tomita finds solutions using elliptical curves.
See companion python code curves_tomita.py

These notes are very messy scribbles as I am learning.
I am slowly cleaning them up.
The goal is that these will be more narrative than the companion python code.
And that this will contain examples from the SageMath envitonment.
Because the expressions are slightly different in sage and in python.

This is a Gemini explanation (But I doubt that it is True):
1. The Weierstrass Transformation
For the surface $x^4 + y^4 + z^4 = 1$, we can define the parameter $m$ as you did. For a specific choice of signs and permutations, the coordinates $(X, Y)$ on the elliptic curve are typically defined as follows:$$X = \frac{2(1 + z^2)}{x^2 + y^2}$$$$Y = \frac{4z(x^2 - y^2)}{(x^2 + y^2)^2}$$
This maps the points on your quartic surface to a Weierstrass equation of the form:$$Y^2 = X^3 - 4X$$(Note: The specific coefficients $A$ and $B$ vary depending on which of the 12 rational "mappings" $m$ you choose.)

This outlines 5 curves and their derivation:
https://math.stackexchange.com/questions/4852495/still-more-elliptic-curves-for-a4b4c4-d4

This outlines Bremner's and Durman's curves:
https://math.stackexchange.com/questions/4853727/on-counterexamples-to-eulers-conjecture-using-bremners-and-durmans-elliptic-c

These give tables by Tomita of curves and solutions:
http://www.maroon.dti.ne.jp/fermat/dioph225e.html
http://www.maroon.dti.ne.jp/fermat/grouplaw1e.html

This gives Tomita's most explicit explanation of his method;
He is following Elkies' explanation.
http://www.maroon.dti.ne.jp/fermat/dioph4e.html

Here is the general outline from that:
2.How to solve the equation A4+B4+C4=D4
    
　　　　It is same as to find the rational solutions of r4+s4+t4=1
        m,n:integer
      　x,y,r,s,t:rational
      　r=x+y,s=x-y

        (2m2+n2)y2=-(6m2-8mn+3n2)x2-2(2m2-n2)x-2mn................ (1) 

        (2m2+n2)t2=4(2m2-n2)x2+8mnx+(n2-2m2)...................... (2)
        
　　　　1. Choose (m,n), and find the rational solution (x,y) of  (1).
　　　　2. Parametrize the rational solution (x,y).
　　　　3. Substitute the rational solution x to (2). (It becomes the elliptic curve.).

　　　　4. Find the rational solutions for elliptic curve.
　　　　5. Convert the rational solutions, and obtain (r,s,t).


Here is how he applies that method to the Frye solution 
(with questions and Gemini explanations in parentheses):
3.Search results

1. Find  958004+2175194+4145604=4224814

　　Frye has already found this solution,but his method is different from Elkies.
    So,I will find a solution  by the method of Elkies theoretically.
　　Elkies pointed out that it was found with (m,n)=(20,-9)
　　Substitute (m,n)=(20,-9) to (1),(2)
　　

　　　　　881y2 = -4083x2-1438x+360........................... (3)

　　　　　881t2 = 2876x2-1440x-719............................ (4)
    Find solution for (3),and obtain (x,y)=(49/318,23/106), and parametrize it
(How find solution x, y?)
sage: R.<x, y, z> = QQ[]
sage: f = 881*y^2 + 4083*x^2 + 1438*x*z - 360*z^2
sage: C = Conic(f)
sage: print(f"Conic: {C}")
Conic: Projective Conic Curve over Rational Field defined by 4083*x^2 + 881*y^2 + 1438*x*z - 360*z^2
sage: pts = C.rational_points(bound=400)
sage: pts
[(1/354 : -75/118 : 1),
 (1/354 : 75/118 : 1),
 (49/318 : -23/106 : 1),
 (49/318 : 23/106 : 1)]
sage: P = C([QQ(49)/318, QQ(23)/106, 1])
sage: P in C
True
phi = C.parametrization(P)
sage: phi
(Scheme morphism:
   From: Projective Space of dimension 1 over Rational Field
   To:   Projective Conic Curve over Rational Field defined by 4083*x^2 + 881*y^2 + 1438*x*z - 360*z^2
   Defn: Defined on coordinates by sending (x : y) to
         (342396*x^2 - 243156*x*y + 43169*y^2 : -491280*x^2 + 345638*x*y - 60789*y^2 : 2229134*x^2 - 1580514*x*y + 280158*y^2),
 Scheme morphism:
   From: Projective Conic Curve over Rational Field defined by 4083*x^2 + 881*y^2 + 1438*x*z - 360*z^2
   To:   Projective Space of dimension 1 over Rational Field
   Defn: Defined on coordinates by sending (x : y : z) to
         (-1208397/3973682*x + 172819/3973682*y + 111849/1986841*z : -1500375107/1750406921*x + 245640/1986841*y + 278152758/1750406921*z))

But this Projective result differs from Tomia's Affine result:
　　      x=1/318(43169*k^2-657351-121578k)/(881k^2+4083)....... (5)

In Sage, the "From" part shows the mapping. It sends a point $(x : y)$ on a line to a point $(X : Y : Z)$ on your curve.To convert Sage's result to the author's "x" (which is $X/Z$):Set Sage's $y = 1$.Set Sage's $x = k$.Calculate $X/Z$.Using your Sage output:Sage's X: $342396k^2 - 243156k + 43169$Sage's Z: $2229134k^2 - 1580514k + 280158$If you divide these, you get a function for $x$ in terms of $k$.

sage: R_k.<k> = QQ[]
sage: x0 = QQ(49)/318
sage: y0 = QQ(-23)/106
sage: var('x_var')
x_var
sage: line = k*(x_var - x0) + y0
sage: eq = 881*(line)^2 + 4083*x_var^2 + 1438*x_var - 360
sage: roots = solve(eq == 0, x_var)
sage: x_k = roots[1].rhs().simplify_full()
sage: x_k
1/318*(43169*k^2 + 121578*k - 657351)/(881*k^2 + 4083)

          y=-1/106(20263*k2-93909+285806k)/(881k2+4083)....... (6)
(Now, $x$ and $y$ are no longer fixed numbers; they are formulas that represent every possible rational point on that conic.)

   Substitute (5) to (4) 
        
         Y2=19435071440k4-5351620404k3+130338882000k2-194951575764k-357457601448.... (7)
         Y=t(140079k2+649197)
(This new equation is a quartic in $k$: $Y^2 = \text{Polynomial of degree 4 in } k$.A quartic equation of the form $Y^2 = f(k)$ is birationally equivalent to an Elliptic Curve in Weierstrass form.)

    Find rational solution for (7),and obtain k=-59/81,Y=2444043484/6561
("Find rational solutions" here means using Elliptic Curve Addition (chord-and-tangent).)

    Substitute k=-59/81 to (5),(6),then x=-159380/422481,y=85060/140827
(The author finds $k = -59/81$. This value of $k$ isn't just a random guess; it is usually found by looking for points of small height or by transforming Equation (7) into a standard Weierstrass $Y^2 = X^3 + AX + B$ and using software (like Magma or Sage) to find generators of the group of rational point)

    r=x+y=-159380/422481+85060/140827=95800/422481
    
    s=x-y=-159380/422481-85060/140827=-414560/422481

    t(140079k2+649197)=2444043484/6561
   
    t=217519/422481

    (95800/422481)4+(-414560/422481)4+(217519/422481)4=1　

    Consequently
　　958004+2175194+4145604=4224814　

More explanation from Gemini on how to find k= -59/81:

1. The Mordell-Weil Group and GeneratorsRational points on an elliptic curve $E$ form an abelian group. This means that if you have two points $P$ and $Q$ on the curve, you can "add" them (using the chord-and-tangent method) to get a third point $R$.The group is structured as:$$E(\mathbb{Q}) \cong T \oplus \mathbb{Z}^r$$$T$ (Torsion): A finite set of points that eventually loop back to the identity (e.g., $P + P + P = \mathcal{O}$).$\mathbb{Z}^r$ (Free part): This is where the infinite solutions live. The Rank ($r$) tells you how many "independent directions" you can move in.Generators: These are the base points $\{P_1, P_2, \dots, P_r\}$. Every possible rational solution $k$ can be expressed as a linear combination: $k(P) = n_1P_1 + n_2P_2 + \dots + \text{Torsion}$.

# Define the curve from the author's Equation (7)
# Note: Equation (7) must be converted to Weierstrass form Y^2 = X^3 + AX + B first
E = EllipticCurve([A, B]) 

# Find the rank and generators
rank = E.rank()
gens = E.gens() 
print(f"Rank: {rank}, Generators: {gens}")

2-DescentImagine the group $E(\mathbb{Q})$. A 2-descent looks at the quotient group $E(\mathbb{Q})/2E(\mathbb{Q})$.It effectively searches for points $P$ such that $P = 2Q$ might be a point on the curve.It maps the curve into several "auxiliary" curves (called homogeneous spaces). These spaces are easier to search for small points.If you find a point on a homogeneous space, you can "map it back" to the original curve $E$. This point is often a generator or a multiple of one.

4-DescentSometimes 2-descent isn't enough because the "gap" between the auxiliary curve and the original curve is still too large.4-descent is a "descent on a descent." It breaks the problem down further into even more auxiliary curves that are geometrically simpler (often represented as intersections of quadrics in $\mathbb{P}^3$).It is computationally expensive but can find generators that are otherwise "invisible" to 2-descent.

# Perform 2-descent to find the rank bounds
rank_info = E.simon_two_descent()

# Search for points using descent methods
points = E.rational_points(bound=10**10) # Simple search
# Or more advanced:
from sage.schemes.elliptic_curves.gp_wrapper import gp
# Use specialized algorithms for higher descent

3. How the Author found $k = -59/81$The author likely followed this path:Converted Equation (7) (the quartic $Y^2 = f(k)$) into the standard Weierstrass form $y^2 = x^3 + Ax + B$.Ran a Descent: Used a 2-descent (likely via a tool like mwrank or Sage's EllipticCurve) to prove the rank was at least 1.Identified the Generator: The descent identified a generator point $P$.Mapped Back: Transformed the $X$ coordinate of the generator $P$ back through the birational maps to the variable $k$.Verified: Plugged $k = -59/81$ into the original equations to get the Frye solution.

Why are there 12 solutions in your table?If the rank of the Frye curve ($u = -9/20$) is $r=1$, the points in your table are likely just the multiples of the generator: $P, 2P, 3P, \dots, 12P$. Because the "height" (size of the numbers) grows quadratically as you take multiples ($nh \approx n^2 h$), you eventually hit the $10^{27}$ limit.

To transform a quartic of the form $Y^2 = f(k)$ into a Weierstrass form $y^2 = x^3 + Ax + B$, SageMath uses a birational transformation. This is essential because the standard algorithms for Descent and Generators only work on the Weierstrass form.

# The coefficients from the author's Equation (7)
# Y^2 = f(k) = ak^4 + bk^3 + ck^2 + dk + e
a = 19435071440
b = -5351620404
c = 130338882000
d = -194951575764
e = -357457601448

# Define the polynomial
R.<k> = QQ[]
f = a*k^4 + b*k^3 + c*k^2 + d*k + e

# Transform to Weierstrass form
# We need one rational point to start the map. 
# The author found k = -59/81. Let's use that.
P_start = [-59/81, 2444043484/6561]
E_trans = EllipticCurve_from_quartic(f, P_start)

print(E_trans)
# This will output something like: Elliptic Curve defined by y^2 = x^3 + Ax + B

2. Finding the GeneratorOnce you have E_trans, you can ask Sage to find the points. The point $k = -59/81$ likely corresponds to a generator $P$ or a small multiple of it.

# Find the generators of the group
G = E_trans.gens()
print(f"Generator Point: {G[0]}")

# If you want to see the 12 solutions in your table, 
# you can iterate through multiples of the generator:
for i in range(1, 13):
    print(f"Point {i}: {i * G[0]}")

Gemini said
To transform a quartic of the form Y 
2
 =f(k) into a Weierstrass form y 
2
 =x 
3
 +Ax+B, SageMath uses a birational transformation. This is essential because the standard algorithms for Descent and Generators only work on the Weierstrass form.

1. SageMath Code for Transformation
You can use the EllipticCurve_from_quartic function. This handles the messy algebra of mapping the coefficients of your Equation (7) to the A and B coefficients.

Python
# The coefficients from the author's Equation (7)
# Y^2 = f(k) = ak^4 + bk^3 + ck^2 + dk + e
a = 19435071440
b = -5351620404
c = 130338882000
d = -194951575764
e = -357457601448

# Define the polynomial
R.<k> = QQ[]
f = a*k^4 + b*k^3 + c*k^2 + d*k + e

# Transform to Weierstrass form
# We need one rational point to start the map. 
# The author found k = -59/81. Let's use that.
P_start = [-59/81, 2444043484/6561]
E_trans = EllipticCurve_from_quartic(f, P_start)

print(E_trans)
# This will output something like: Elliptic Curve defined by y^2 = x^3 + Ax + B
2. Finding the Generator
Once you have E_trans, you can ask Sage to find the points. The point k=−59/81 likely corresponds to a generator P or a small multiple of it.

Python
# Find the generators of the group
G = E_trans.gens()
print(f"Generator Point: {G[0]}")

# If you want to see the 12 solutions in your table, 
# you can iterate through multiples of the generator:
for i in range(1, 13):
    print(f"Point {i}: {i * G[0]}")
3. Understanding the Descent Process
When Sage executes E.gens(), it often runs a 2-Descent behind the scenes. Here is what is actually happening algebraically:

2-Descent (The "Half-Point" Search)
The goal is to find the rank. The algorithm looks for points P that are "double" another point (P=2Q).

It checks if the equation y 
2
 =f(x) can be factored into auxiliary equations called 2-covers.

If a 2-cover has a rational point, it contributes to the rank.

Sage Notation: E.rank_bounds() or E.two_descent().

4-Descent (The "Quarter-Point" Search)
If the 2-descent is inconclusive (it gives a range, like "rank is between 0 and 2"), we use 4-descent.

It lifts the 2-covers to 4-covers (intersections of two quadrics in 3D space).

It is much harder to find points on a 4D quadric intersection, but if you find one, it almost certainly leads to a generator of the original curve.

Sage Notation: E.four_descent().

4. What it means for your 12 solutionsThe fact that you have 12 solutions for $u = -9/20$ suggests that:The Frye curve has Rank 1 (one infinite generator $P$).The points $(a, b, c, d)$ in your table correspond to the points $nP$ for $n = 1, 2, \dots, 12$ on that specific elliptic curve.Because the Canonical Height $h(nP) = n^2 h(P)$, the number of digits in your $d$ value grows quadratically. If $P$ has a $d$ around $10^2$, then $12P$ would have a $d$ around $10^{2^2 \times 12^2} \dots$ effectively reaching your $10^{27}$ limit very quickly.

I gave Gemini the 12 solutions on curve C1 from my piezas_in_pairs.py

1. Are these multiples of a single generator?Almost certainly. When an elliptic curve has Rank 1, every rational point is of the form $nP$.Height Growth: On an elliptic curve, the "size" of the coordinates (the logarithmic height) grows quadratically: $Height(nP) \approx n^2 \times Height(P)$.The Bound: Your largest $d$ is approximately $10^{26}$. If the generator $P$ is the Frye solution ($d \approx 10^5$), then $7P$ or $8P$ would already be pushing the limits of $10^{27}$.Tomita's Solutions: Notice that many of the later solutions (like #80, #84, #87) were found in 2024. These represent the result of massive computational searches (often using 4-descent or high-performance sieving) to find the "next" multiples in the sequence that were previously too large to identify.

2. Why do some share other rationals?Look at your first few entries:Solution 1 (Frye): Shares $u \in \{-1041/320, 1000/47\}$.Solution 17 (Tomita): Also shares $u \in \{-1041/320, 1000/47\}$.This means Solution 1 and Solution 17 are likely related by a very simple operation in the group law. If Solution 1 is the generator $P$, Solution 17 might be $-P$ or $2P$. Because the 12 rationals are computed by permuting $a, b, c$, a point and its "negative" in the group law will often share a large subset of these $u$-values, as the coordinates only differ by signs or specific rational shifts.

Gemini said
This data is a goldmine for understanding how the group law of elliptic curves works in practice. By looking at the d values and the sets of rationals (the u-values), we can map out exactly how these points relate to each other on the "Frye Curve" (u=−9/20).

1. Are these multiples of a single generator?
Almost certainly. When an elliptic curve has Rank 1, every rational point is of the form nP.

Height Growth: On an elliptic curve, the "size" of the coordinates (the logarithmic height) grows quadratically: Height(nP)≈n 
2
 ×Height(P).

The Bound: Your largest d is approximately 10 
26
 . If the generator P is the Frye solution (d≈10 
5
 ), then 7P or 8P would already be pushing the limits of 10 
27
 .

Tomita's Solutions: Notice that many of the later solutions (like #80, #84, #87) were found in 2024. These represent the result of massive computational searches (often using 4-descent or high-performance sieving) to find the "next" multiples in the sequence that were previously too large to identify.

2. Why do some share other rationals?
Look at your first few entries:

Solution 1 (Frye): Shares u∈{−1041/320,1000/47}.

Solution 17 (Tomita): Also shares u∈{−1041/320,1000/47}.

This means Solution 1 and Solution 17 are likely related by a very simple operation in the group law. If Solution 1 is the generator P, Solution 17 might be −P or 2P. Because the 12 rationals are computed by permuting a,b,c, a point and its "negative" in the group law will often share a large subset of these u-values, as the coordinates only differ by signs or specific rational shifts.

3. Verification via SageMath
To prove these are all on the same "line" (the same cyclic group), we use the transformation we discussed. We can take the coordinates from Solution 1 and Solution 2, map them to (X,Y) coordinates on the Weierstrass curve, and check if P2 =nP1

# SageMath logic to check the relationship
# Define the curve based on the -9/20 parameter
E = EllipticCurve([0, A, 0, B, C]) # Using the transformed Weierstrass coefficients

# Map Frye's solution (Sol 1) to a point on the curve
P1 = E(x1, y1) 

# Map MacLeod's solution (Sol 2) to a point on the curve
P2 = E(x2, y2)

# Check the relationship
for n in range(-10, 10):
    if n*P1 == P2:
        print(f"Solution 2 is the {n}-th multiple of Solution 1")

This post by Tito Piezas explores the structure of the curves.
https://math.stackexchange.com/questions/4852495/still-more-elliptic-curves-for-a4b4c4-d4

For example it gives the 12 v paired with u=-9/20 to find the 12 solutions.

Dan Fulea also talks about hyperelliptical curves. Gemini explains:
Feature,Elliptic Curve,Hyperelliptic Curve (Genus 2)
Equation,y2=cubic,y2=quintic/sextic
Number of Points,Can be infinite (if Rank > 0),Always finite (Faltings' Theorem)
Sage Class,EllipticCurve(f),HyperellipticCurve(f)
Group Law,Points on the curve form a group,Points on the Jacobian form a group

In this essay, Tomita analyzes Elkie's (8,-5)
http://www.maroon.dti.ne.jp/fermat/grouplawe.html

This is Gemini's answer to some of my questions about the essay:

Gemini said
The author is using a sophisticated toolset that bridges abstract geometry and high-performance computation. Here is the breakdown of those technical components.

1. Transforming to Minimal Weierstrass Form
To move from the quartic V2 =f(U) to the Weierstrass form y2 =x3 +Ax+B, the author uses a birational map.
The Process: This involves picking a rational point on the quartic and moving it to the "point at infinity" on the elliptic curve.
Minimal Form: A "Minimal" Weierstrass form is one where the coefficients A and B are integers and are as small as possible (specifically, the discriminant Δ has no factors of p^12
  that can be removed by a change of variables). This makes the computation of rank and points significantly faster.

2. Sage rank() vs. mwrank
You likely found that Sage's default E.rank() hung because it often tries to perform a full descent using general-purpose algorithms.

What is mwrank? It is a specialized C++ program written by John Cremona. It is highly optimized for 2-descent.

Is it more feasible? Yes. For curves with large coefficients or higher ranks (like this Rank 3 curve), mwrank is the industry standard. It uses clever sieving and modular arithmetic to search the "homogeneous spaces" created during descent much faster than general Sage methods.

How it works: It breaks the curve into several simpler "covers." If it finds a point on any of these covers, it maps it back to the main curve.

3. Understanding the NotationThe Colons (Projective Coordinates)The expressions like [40984562244 : 7472357686746 : 1225043] are Projective Coordinates $[X : Y : Z]$.In projective space, we avoid fractions by adding a third coordinate ($Z$).To get the standard $(x, y)$ coordinates you recognize, you divide by $Z$ (or a power of $Z$):$$x = \frac{X}{Z}, \quad y = \frac{Y}{Z}$$Example: $x = 40984562244 / 1225043 \approx 33455.6 \dots$

The colons (e.g., [15651603708098045508 : 5912675769930792830522 : 123885826842051]) are just a way of writing the fraction $X/Z^2$ and $Y/Z^3$ without actually performing the division. Because the height of these points is so large ($26+$), the number of digits in the $X, Y, Z$ coordinates will naturally be huge.

The Canonical Height (or Weil Height) is a measure of the "complexity" or "size" of a point.It is roughly proportional to the number of digits in the numerator and denominator.Why it matters: In the search for points, $mwrank$ searches in order of increasing height. A height of $28.14$ (Generator 3) means the numbers involved are massive (roughly $e^{28}$), which is why that point is so much harder to find than Generator 2 (height $11.15$).

The Regulator is a value that measures the "density" or "spread" of the generators.

Mathematically, it is the determinant of a matrix built from the heights of the generators.

If the regulator is small, the points are "close" together in terms of height. If it is large, the points are very spread out. It is a key component in the Birch and Swinnerton-Dyer (BSD) Conjecture, which relates the number of points to the analytical properties of the curve.

4. Why two different descriptions for the points?You noticed the author first lists the points with colons, and then as fractions.The mwrank Output (Colons): This is the raw output of the software. mwrank works in projective space $[X:Y:Z]$ because it is computationally more stable to use large integers than to deal with floating-point or rational fraction overhead during the descent phase.The Rational Points (Fractions): The author later converts them to the form $P = [x, y]$ to make them usable in the formulas for $k, x, y, r, s$.Note: In the second list, there is a typo in the author's text: p1=[[383033292/11449.... The denominator $11449$ is actually the square root of the $Z$ part if they were mapping back from a different model, or simply the simplified fraction form of $X/Z$.

# Define the curve from the essay
# y^2 = x^3 - x^2 + 2815805388*x - 94443526967868
a1, a2, a3, a4, a6 = 0, -1, 0, 2815805388, -94443526967868
E = EllipticCurve([a1, a2, a3, a4, a6])

# Use mwrank via Sage's interface
from sage.libs.cremona.all import mwrank_curve

# This initializes the mwrank engine for this specific curve
m = mwrank_curve([a1, a2, a3, a4, a6])

# Get the rank and generators
print(f"Rank: {m.rank()}")
print(f"Generators: {m.gens()}")

But that gave an error failing to import mwrank_curve.
Instead:
# Define the curve coefficients
# a1, a2, a3, a4, a6
coeffs = [0, -1, 0, 2815805388, -94443526967868]

# Call mwrank directly
# verbose=False keeps the output clean, but True shows the descent progress
results = mwrank(coeffs)
print(results)

I got slightly different results than Tomita 
(He had Generator P3 = E([287246048321418303786, 267958777554236828867184, 331273336732057])  Height $\approx 28.15$  )
I get
Generator 1 is [40984562244:7472357686746:1225043]; height 13.048394176944
Generator 2 is [7682852:2843134723:64]; height 11.157481098979
Generator 3 is [1053454466437509794928036:258481066147406623852872939:23635653014313151552]; height 33.937039835591
Regulator = 1267.61149551
The rank and full Mordell-Weil basis have been determined unconditionally.
 (0.241227 seconds)

How to use the curve:
E = EllipticCurve(coeffs)
Elliptic Curve defined by y^2 = x^3 - x^2 + 2815805388*x - 94443526967868 over Rational Field

# Manually define the generators from the author's mwrank output
# Sage handles the [X : Y : Z] projective notation automatically!
P1 = E([40984562244, 7472357686746, 1225043])
(383033292/11449 : 7472357686746/1225043 : 1)
sage: P1 in E
True
same for P2 and his Generator 3 and mine.

# Example: Finding a new solution by adding two generators
new_point = P1 + P2
(42396 : -10056690 : 1)
print(f"New Point x-coordinate: {new_point[0]}")
2*new_point
(100135114/1225 : 1119621383712/42875 : 1)

If mwrank fails to finish, use options
        -h       help            prints this info and quits
 |         -v n     verbosity       sets verbosity to n (default=1)
 |         -o       PARI/GP output  turns ON extra PARI/GP short output (default: OFF)
 |         -p n     precision       sets precision to n decimals (default=15)
 |         -b n     quartic bound   bound on quartic point search (default=10)
 |         -x n     n aux           number of aux primes used for sieving (default=6)
 |         -l       list            turns ON listing of points (default: ON, unless v=0)
 |         -s       selmer_only     if set, computes Selmer rank only (default: not set)
 |         -d       skip_2nd_descent        if set, skips the second descent for curves with 2-torsion (default: not set)
 |         -S n     sat_bd          upper bound on saturation primes (default=100, -1 for automatic)
For example,
M = Mwrank('-b 5 -p 20')
coeffs = [0, -1, 0, 2815805388, -94443526967868]
output = M(coeffs)
print(output)
This time got another new 3rd Generator
Generator 3 is [15651603708098045508:5912675769930792830522:123885826842051]; height 26.519
The difference is that the generators are just different basis for same group.
The group is the same because the regulator is unchanged.

Convert my Generator 3 to a point:
sage: X, Y, Z = 15651603708098045508, 5912675769930792830522, 123885826842051
sage: E([X, Y, Z])
(313967697901708/2485122201 : 5912675769930792830522/123885826842051 : 1)

2. Why mwrank gives different resultsThe mwrank algorithm (specifically the 2-descent) works by searching through "quartic covers."Search Order: It searches for points based on their height (complexity).The "Quartic Bound" (-b): When you change the search bound or the precision, mwrank might find points in a different order.Saturation: Once mwrank finds enough independent points to satisfy the rank (3 in this case), it stops. If your search parameters or the internal state of the library shifted, it might "grab" a different point that is just as valid as the author's $P_3$.

my_gens = [P1, P2, P3_yours]
# Create a list of the author's points
author_gens = [P1, P2, P3_author]

# Check the regulator (Sage command)
E = EllipticCurve([0, -1, 0, 2815805388, -94443526967868])
print(E.regulator_of_points(my_gens))
print(E.regulator_of_points(author_gens))

Since my last Generator 3 has smaller height than the others, it might be able to find
a new solution, but since it is higher than 1 and 2, it would probably be > 1e27.

This runs the mwrank logic and returns the generators as points:
E.gens()
[(383033292/11449 : 7472357686746/1225043 : 1),
 (367079304556757337/8235912749584 : 258481066147406623852872939/23635653014313151552 : 1),
 (1920713/16 : 2843134723/64 : 1)]

M = Mwrank('-b 20 -p 30')
output = M(coeffs)
E.gens()
same result.

The table post:
https://math.stackexchange.com/questions/1853223/distribution-of-primitive-pythagorean-triples-ppt-and-of-solutions-of-a4b4/4857107#4857107

lists the 5 most prolific curves
-5/8 𝐶0:𝑋3−𝑋2+2815805388𝑋−94443526967868=𝑌2
analyzed above.

-9/20 𝐶1:𝑋3+2265722465761𝑋−3154189403034549278=𝑌2
coeffs = [0, -1, 0, 2265722465761, -3154189403034549278]
results = mwrank(coeffs)
I gave up waiting
Review posts for what Tomita did when couldn't get generators.
Or ask Gemini how to coax it.

-41/36 𝐶5:𝑋3+2639323244332897𝑋−20156152630838819347102=𝑌2
coeffs = [0, -1, 0, 2639323244332897, -20156152630838819347102]
I gave up waiting

-5/44 𝐶7:𝑋3−𝑋2+349942184229228𝑋−11167797929528591502588=𝑌2
coeffs = [0, -1, 0, 349942184229228, -11167797929528591502588]

-125/92 𝐶9:𝑋3−𝑋2+9243195710310751148𝑋−761969307339454319105751548=𝑌2
coeffs = [0, -1, 0, 9243195710310751148, -761969307339454319105751548]


Would you like me to write the formula to transform a point $(x, y)$ from this curve back into the $A, B, C, D$ solution using the author's $U \to k \to A, B, C, D$ steps?

The author uses a chain of variables. To go backward, we follow this sequence:
Point $(x, y)$ on (4) $\to$ Coordinate $U$ on (3) $\to$ Parameter $k$ $\to$ Ratios $x, y$ $\to$ Final $r, s, t$.

def point_to_solution(P):
    # 1. Map Weierstrass x to Quartic U
    # This mapping depends on the specific birational transformation used.
    # For the (8, -5) case, U is a rational function of the Weierstrass x.
    # Note: If you used 'EllipticCurve_from_quartic', Sage stores the map!
    
    # Assuming we have the U value (the author's quartic variable):
    U = P[0] # Simplification: In many mappings, U is linearly related to x
    
    # 2. Author's Step: k = U*7 - 2
    k = U * 7 - 2
    
    # 3. Author's formulas for x and y (the internal ratios)
    den_base = 14 * (17*k^2 + 779)
    x_ratio = (51*k^2 - 34*k - 5221) / den_base
    
    den_base_y = 42 * (17*k^2 + 779)
    y_ratio = (17*k^2 + 7558*k - 779) / den_base_y
    
    # 4. Convert to r, s, t (A/D, B/D, C/D)
    r = x_ratio + y_ratio
    s = x_ratio - y_ratio
    
    # 5. Solve for t using r^4 + s^4 + t^4 = 1
    t4 = 1 - r^4 - s^4
    t = t4^(1/4) # This should be rational if the point is correct
    
    # 6. Integerize to (A, B, C, D)
    lcm_val = r.denominator().lcm(s.denominator()).lcm(t.denominator())
    A = abs(r * lcm_val)
    B = abs(s * lcm_val)
    C = abs(t * lcm_val)
    D = abs(lcm_val)
    
    return sorted([int(A), int(B), int(C)]), int(D)

# Example usage with your P_combined:
# sol_abc, sol_d = point_to_solution(P_combined)
# print(f"{sol_abc[0]}^4 + {sol_abc[1]}^4 + {sol_abc[2]}^4 = {sol_d}^4")

When the author says "Find rational solution $U$ for (3)," he is identifying $U$ as the link between the geometry of the surface and the arithmetic of the curve.The quartic (3) is essentially a "double cover" of the $U$-line.The Weierstrass curve (4) is just a mathematically "standard" version of that quartic.

3. Important Caveat on the $x \to U$ MapThe author omitted the explicit map from his Weierstrass $x$ to his quartic $U$. However, when you use Sage's EllipticCurve_from_quartic function, Sage creates an object that knows the map:

# If you created E this way:
E, f_map, inverse_map = EllipticCurve_from_quartic(quartic_poly, return_maps=True)

# You can get U directly from a point P on E:
U = inverse_map(P)[0]



Experiments on tuvd_to_u_list
m = (a^2+2ab+b^2-d^2-c^2)/(a^2+ab+b^2-da-db)

#1 tuvw_to_u_list((11_975, 51_820, 102_481, 320_000))
95800 414560 217519 422481

(52463, 660460, 660706)
(-4209, 3500, 59914800)
(3521543, 9580960, 60160)

(52463, 660460, 660706) dupe
(-1000, -47, 335534674)
(-71490240, 101943281, 2082)

(-6899820729, 369596780, 18)
(-4209, 3500, 59914800) dupe
(-1041, 320, 142980480)

(-6899820729, 369596780, 18) dupe
(-1000, -47, 335534674) dupe
(30080, 6007, 7043086)

(-9, 20, 13799641458)
(-167767337, 43538900, 2000)
(3521543, 9580960, 60160) dupe

(-9, 20, 13799641458) dupe
(-29957400, -6538471, 8418)
(-71490240, 101943281, 2082) dupe

(330353, 48940, 104926)
(-167767337, 43538900, 2000) dupe
(-1041, 320, 142980480) dupe

(330353, 48940, 104926) dupe
(-29957400, -6538471, 8418) dupe
(30080, 6007, 7043086) dupe
(24, 12)

#2 tuvw_to_u_list((173_800, 345_953, 1_069_568, 1_743_433))
1390400 2767624 673865 2813001

(130344703, 367240100, 68450)
(-41952, 33865, 269630226)
(105512, 1038385, 18933010)

(130344703, 367240100, 68450) dupe
(-1314533417, 302546960, 11456)
(-9, 20, 606928892880)

(-1425, 412, 4540577346)
(-41952, 33865, 269630226) dupe
(-303464446440, 19949196449, 18)

(-1425, 412, 4540577346) dupe
(-1314533417, 302546960, 11456) dupe
(9466505, 1411828, 211024)

(-2270288673, 3375081700, 2850)
(-5728, -215, 2629066834)
(105512, 1038385, 18933010) dupe

(-2270288673, 3375081700, 2850) dupe
(-134815113, -29587760, 83904)
(-9, 20, 606928892880) dupe

(34225, 6692, 260689406)
(-5728, -215, 2629066834) dupe
(-303464446440, 19949196449, 18) dupe

(34225, 6692, 260689406) dupe
(-134815113, -29587760, 83904) dupe
(9466505, 1411828, 211024) dupe
(24, 12)

#3 tuvw_to_u_list((688_485, 1_041_526, 3_500_953, 5_206_528))
5507880 8332208 1705575 8707481

(103611, 244444, 1088869178)
(-7421336185, 8395768548, 12560)
(-3333, 107368, 1620125984)

(103611, 244444, 1088869178) dupe
(-35063912120, 15270810801, 3730)
(-30768, 57253, 2023033754)

(-29, 12, 2439737337738)
(-7421336185, 8395768548, 12560) dupe
(-1011516877, 6860856, 61536)

(-29, 12, 2439737337738) dupe
(-35063912120, 15270810801, 3730) dupe
(-810062992, -128701299, 6666)

(-1219868668869, 1352800221524, 58)
(-1865, -132, 70127824240)
(-3333, 107368, 1620125984) dupe

(-1219868668869, 1352800221524, 58) dupe
(-6280, -1359, 14842672370)
(-30768, 57253, 2023033754) dupe

(544434589, 121333308, 207222)
(-1865, -132, 70127824240) dupe
(-1011516877, 6860856, 61536) dupe

(544434589, 121333308, 207222) dupe
(-6280, -1359, 14842672370) dupe
(-810062992, -128701299, 6666) dupe
(24, 12)

#4 tuvw_to_u_list((733_750, 1_411_130, 1_957_457, 10_240_000))
5870000 11289040 8282543 12197457

(20632147, 117135680, 3734666)
(-84237, 359800, 901336800)
(8685847, 22963880, 23004320)

(20632147, 117135680, 3734666) dupe
(-400, 37, 676000044154)
(-35798568240, 28609248113, 4866)

(-1010819791893, 158785763120, 186)
(-84237, 359800, 901336800) dupe
(-2433, 920, 71597136480)

(-1010819791893, 158785763120, 186) dupe
(-400, 37, 676000044154) dupe
(11502160, 2925527, 17371694)

(-93, 80, 2021639583786)
(-338000022077, 104832225800, 800)
(8685847, 22963880, 23004320) dupe

(-93, 80, 2021639583786) dupe
(-450668400, -124346123, 168474)
(-35798568240, 28609248113, 4866) dupe

(1867333, 457280, 41264294)
(-338000022077, 104832225800, 800) dupe
(-2433, 920, 71597136480) dupe

(1867333, 457280, 41264294) dupe
(-450668400, -124346123, 168474) dupe
(11502160, 2925527, 17371694) dupe
(24, 12)
(-93, 80, 2021639583786)
(-400, 37, 676000044154)

#5 tuvw_to_u_list((1_569_025, 1_771_715, 5_761_993, 10_241_024))
12552200 14173720 4479031 16003017

(7919435, 17426416, 55321690)
(-511289, 1551044, 326484752)
(-7516797, 66665212, 8744304)

(7919435, 17426416, 55321690) dupe
(-974663766264, 623365202417, 402)
(-136, 133, 2350494728938)

(-1005, 568, 272168104170)
(-511289, 1551044, 326484752) dupe
(-1175247364469, 8556077924, 272)

(-1005, 568, 272168104170) dupe
(-974663766264, 623365202417, 402) dupe
(-4372152, -935219, 15033594)

(-136084052085, 102731297224, 2010)
(-201, -4, 1949327532528)
(-7516797, 66665212, 8744304) dupe

(-136084052085, 102731297224, 2010) dupe
(-163242376, -37855247, 1022578)
(-136, 133, 2350494728938) dupe

(27660845, 6860848, 15838870)
(-201, -4, 1949327532528) dupe
(-1175247364469, 8556077924, 272) dupe

(27660845, 6860848, 15838870) dupe
(-163242376, -37855247, 1022578) dupe
(-4372152, -935219, 15033594) dupe
(24, 12)

#6 tuvw_to_u_list((455_355, 878_575, 74_752, 16_355_761))
3642840 7028600 16281009 16430513

(-956745, 598856, 440193970)
(1722153304, 14758970841, 45058)
(79416, 247889, 3275232082)

(-956745, 598856, 440193970) dupe
(-282464081, 20137068, 565168)
(-7683273, 1866844, 25721456)

(-12185, -432, 42968374290)
(1722153304, 14758970841, 45058) dupe
(-12860728, 22910913, 15366546)

(-12185, -432, 42968374290) dupe
(-282464081, 20137068, 565168) dupe
(1637616041, 289087452, 158832)

(-21484187145, 3803742704, 24370)
(-282584, 755289, 564928162)
(79416, 247889, 3275232082) dupe

(-21484187145, 3803742704, 24370) dupe
(22529, 2988, 3444306608)
(-7683273, 1866844, 25721456) dupe

(-220096985, -45498888, 1913490)
(-282584, 755289, 564928162) dupe
(-12860728, 22910913, 15366546) dupe

(-220096985, -45498888, 1913490) dupe
(22529, 2988, 3444306608) dupe
(1637616041, 289087452, 158832) dupe
(24, 12)

#7 tuvw_to_u_list((335_305, 2_349_595, 2_625_017, 17_990_656))
2682440 18796760 15365639 20615673

(-46747125, 199797104, 4273050)
(-477, 692, 948828100944)
(10012917591, 21562651564, 73392)

(-46747125, 199797104, 4273050) dupe
(-20824, -2003, 29651393722)
(-1688076136, 1490944633, 249058)

(-40143752605885, -2211506178056, 10)
(-477, 692, 948828100944) dupe
(-124529, 68084, 3376152272)

(-40143752605885, -2211506178056, 10) dupe
(-20824, -2003, 29651393722) dupe
(36696, 8687, 20025835182)

(-5, 8, 80287505211770)
(-14825696861, 11130273356, 41648)
(10012917591, 21562651564, 73392) dupe

(-5, 8, 80287505211770) dupe
(-474414050472, -91779396259, 954)
(-1688076136, 1490944633, 249058) dupe

(-2136525, -340912, 93494250)
(-14825696861, 11130273356, 41648) dupe
(-124529, 68084, 3376152272) dupe

(-2136525, -340912, 93494250) dupe
(-474414050472, -91779396259, 954) dupe
(36696, 8687, 20025835182) dupe
(24, 12)

#8 tuvw_to_u_list((270_579, 3_958_640, 1_613_041, 42_697_216))
2164632 31669120 41084175 44310257

(-12065, 12396, 207756778818)
(-3622925, 12225132, 302482864)
(76288487203, 165557531760, 43584)

(-12065, 12396, 207756778818) dupe
(-8583400, -714723, 169117274)
(-594383695904, 307577299845, 3162)

(-1701832006545, -225388398548, 1634)
(-3622925, 12225132, 302482864) dupe
(-1581, 1520, 1188767391808)

(-1701832006545, -225388398548, 1634) dupe
(-8583400, -714723, 169117274) dupe
(21792, 5035, 152576974406)

(-817, 660, 3403664013090)
(-84558637, 193874100, 17166800)
(76288487203, 165557531760, 43584) dupe

(-817, 660, 3403664013090) dupe
(-151241432, -18609075, 7245850)
(-594383695904, 307577299845, 3162) dupe

(-103878389409, -17530520620, 24130)
(-84558637, 193874100, 17166800) dupe
(-1581, 1520, 1188767391808) dupe

(-103878389409, -17530520620, 24130) dupe
(-151241432, -18609075, 7245850) dupe
(21792, 5035, 152576974406) dupe
(24, 12)

#9 tuvw_to_u_list((1_301_137, 5_359_820, 1_389_056, 67_322_041))
10409096 42878560 65932985 68711097

(-68923973, 66997508, 90372026)
(-6552, 93017, 111671615874)
(13531216960, 31827291581, 518090)

(-68923973, 66997508, 90372026) dupe
(-10498601, -138604, 331175888)
(-163997925, 67709216, 26208000)

(-190621456389, -17373563356, 42042)
(-6552, 93017, 111671615874) dupe
(-13104000, 15069437, 327995850)

(-190621456389, -17373563356, 42042) dupe
(-10498601, -138604, 331175888) dupe
(259045, 56768, 27062433920)

(-21021, 9788, 381242912778)
(-165587944, 361205857, 20997202)
(13531216960, 31827291581, 518090) dupe

(-21021, 9788, 381242912778) dupe
(-55835807937, -7916481364, 13104)
(-163997925, 67709216, 26208000) dupe

(-45186013, -9200036, 137847946)
(-165587944, 361205857, 20997202) dupe
(-13104000, 15069437, 327995850) dupe

(-45186013, -9200036, 137847946) dupe
(-55835807937, -7916481364, 13104) dupe
(259045, 56768, 27062433920) dupe
(24, 12)

(-6552, 93017, 111671615874)
(-21021, 9788, 381242912778)

#10 tuvw_to_u_list((4_364_815, 13_270_140, 14_623_232, 102_488_849))
34918520 106161120 87865617 117112081

(-29915585, 638841996, 51215490)
(-15944968, 42458829, 621479994)
(23301824, 52367253, 974680442)

(-29915585, 638841996, 51215490) dupe
(-11846053, -175812, 1872524784)
(-93514757, 75615072, 156123136)

(-431691625, 11587212, 37897650)
(-15944968, 42458829, 621479994) dupe
(-78061568, 40207629, 187029514)

(-431691625, 11587212, 37897650) dupe
(-11846053, -175812, 1872524784) dupe
(487340221, 120067584, 46603648)

(-18948825, 19835764, 863383250)
(-936262392, 509548501, 23692106)
(23301824, 52367253, 974680442) dupe

(-18948825, 19835764, 863383250) dupe
(-310739997, -74371388, 31889936)
(-93514757, 75615072, 156123136) dupe

(-25607745, -5442892, 59831170)
(-936262392, 509548501, 23692106) dupe
(-78061568, 40207629, 187029514) dupe

(-25607745, -5442892, 59831170) dupe
(-310739997, -74371388, 31889936) dupe
(487340221, 120067584, 46603648) dupe
(24, 12)