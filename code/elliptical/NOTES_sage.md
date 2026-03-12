Notes on how to install SageMath and use it.

Install sage as a conda ennvironment.

(base) rogerfrye@MacBookPro code % conda activate sage
(sage) rogerfrye@MacBookPro code % sage
┌────────────────────────────────────────────────────────────────────┐
│ SageMath version 10.7, Release Date: 2025-08-09                    │
│ Using Python 3.12.12. Type "help()" for help.                      │
└────────────────────────────────────────────────────────────────────┘
sage: factor(11_516_209)


Use %hist for command history
Use _oh for output history
Use _i1 to see input 1
Use eval(_i1) to evaluate input 1
Use %macro em 1-2 to create macro em from inputs 1-2.
Use em to excute macro em
Use !ls to execute shell command ls
Use command? for help on command
Use logstart setup to record a session, and later load("setup") to reload it.
Use %time or %timeit to time a computation
Use t = cputime(); ...; cputime(t) to time a block of lines.
time g = gp('1938^99484')
Use %edit to open an editor
Use %pdb to debug
Use %quickref to a man page of these commands
Use E.(tab) after creating object E to see methods
Use help(EllipticCurve) to get help
Use save(x, filename) and load(filename) to save and reload pickled object
	same as cPickle.dumps(x, 2)

There are several ways to construct an elliptic curve
y^2 + a_1 xy + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6 
EllipticCurve([a1,a2,a3,a4,a6])
EllipticCurve([a4,a6])
EllipticCurve(label) # Cremona database of rational curves
EllipticCurve(j=j0) # j-invariant
EllipticCurve(polynomial) # see EllipticCurve_from_Weierstrass_polynomial
EllipticCurve(cubic, point) # (homogeneous polynomial in three variables, rational point

Np(p) The number of points on E modulo p
E1 = EllipticCurve([2265722465761,-3154189403034549278])
E1.conductor()
112191513909404296376
E.Np(97)
116

E1.rank(algorithm="pari")
3
E1.gens(algorithm="pari", pari_effort=15)
[(47971729/49 : 16603172706/343 : 1),
 (1237921 : 1244044242 : 1),
 (98822808873914209/20881117009 : 32155995292015599252526194/3017384051151527 : 1)]

p1 = E1(978559, 0)
p2 = E1(47971729/49, 16603172706/343)
p3 = E1(1237921, 1244044242)
p1+p3
(1684027681/81 : -69276074469938/729 : 1)
2*(p1+p3)
(998224478479/194481 : -1028117289506708080/85766121 : 1)

E.tor Tab. pops up [method for method in dir(E) if 'torsion' in method]
E.torsion_subgroup? for documentation
E.torsion_subgroup?? for code
help(E) for documentation page
search_doc("torsion")	Searches all of Sage's offline documentation for the word "torsion".

See this doc for compile code:
https://doc.sagemath.org/html/en/tutorial/programming.html#section-loadattach

import sage.schemes.elliptic_curves as EC
EC.__path__
_NamespacePath(['/Users/rogerfrye/miniforge/envs/sage/lib/python3.12/site-packages/sage/schemes/elliptic_curves'])

