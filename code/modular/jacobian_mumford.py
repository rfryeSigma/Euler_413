"""
The function u_to_quartic(u), defines a curve of the form
D^2 = f(v), where f(v) is a quartic polynomial. In algebraic geometry, 
a non-singular curve of genus 1 is essentially an elliptic curve without 
a designated origin. The Jacobian of this curve is the unique elliptic 
curve that is birationally equivalent to it and possesses a rational 
point to serve as the identity (O).
The Jacobian is the elliptic curve: Y^2 = X^3 - 27 I X - 27 J
where I and J are the semi-invariants of the quartic.

Mumford representation is usually overkill on the Genus 1 quartic
because we have the standard (x, y) group law.
However, if the coefficients are massive, working with the quartic directly 
(using "Jacobi Quartic" addition formulas) 
can sometimes be more numerically stable than Weierstrass coordinates.

The Fermat Quartic surface can be mapped to the Jacobian of a 
Genus 2 curve y^2 = f(x) where deg(f)=5 or 6. 
This is where Mumford representation u(x), v(x) becomes powerful. 
Adding points in a Genus 2 Jacobian allows you to find solutions 
that aren't just "multiples" of a single u, 
but rather complex combinations that might "jump" between different u.

------------
I studied
Jacobian elliptic Kummer surfaces and special function identities 
    by Elise Griffin and Andreas Malmendier
    Communications in Number Theory and Physics
https://intlpress.com/site/pub/files/_fulltext/journals/cntp/2018/0012/0001/CNTP-2018-0012-0001-a004.pdf

They derive formulas for the construction of all inequivalent Jacobian 
elliptic fibrations on the Kummer surface of two non-isogeneous elliptic 
curves from extremal rational elliptic surfaces by rational base 
transformations and quadratic twists.

-----------
I studied
Kummer Surfaces, Isogenies and Theta Functions
    by Adrian Clingher and Andreas Malmendier
https://arxiv.org/html/2505.13727v2
These transforms are so hard they are being used for cryptography.

TODO
Study:

Automorphisms of Jacobian Kummer surfaces
    by JONG HAE KEUM

Mumford representation.
Jacobi Quartic addition.
Picard number
Picard Group
Picard Lattice on K3.

"""