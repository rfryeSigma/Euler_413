"""
Tito Piezas II identified 5 rank 3 Elliptical Curves
responsible for more than a third of the known solutions.
See
https://math.stackexchange.com/questions/1853223/distribution-of-primitive-pythagorean-triples-ppt-and-of-solutions-of-a4b4
Investigate them.
I have the a_invariants and conductor. No cremona_label.
"""
C0 = EllipticCurve((0, 0, 0, 2815805388, 94443526967868))
C1 = EllipticCurve((0, 0, 0, 2265722465761, 3154189403034549278))
C5 = EllipticCurve((0, 0, 0, 2639323244332897, 20156152630838819347102))
C7 = EllipticCurve((0, 0, 0, 349942184229228, 11167797929528591502588))
C9 = EllipticCurve((0, 0, 0, 9243195710310751148, 761969307339454319105751548)) 

C0.conductor()
110044011515897869948876416912
C1.conductor()
224383027818808592752
C5.conductor()
1480359324892946485840
C7.conductor()
1179615736971507255354138336942050718545426832
C9.conductor()
3174507458503050356219986867886964719461924390925486515376

"""
But LMFDB https://www.lmfdb.org/EllipticCurve/Q/
conductor must be < 300_000 .
"""

"""
What does discriminant() tll me?
$\Delta < 0$: The curve has only one component. 
This is a single infinite curve. 
This happens when the underlying cubic has only one real root and two complex ones.
"""
C0.discriminant()
-5282112552763097757546068011776
C1.discriminant()
-5042318735472245059970154044172491183872
C5.discriminant()
-1352187089796352374424826935630453193995900000000
C7.discriminant()
-56621555374632348256998640173218434490180487936
C9.discriminant()
-50792119336048805699519789886191435511390790254807784246016

"""
When I ask for rank, it churns overnight and doesn't return an answer.
"""
C0.rank()

C1.rank()

C5.rank()

C7.rank()

C9.rank()

"""
What is this stuff that copilot threw out?
C0.torsion_subgroup()
Torsion subgroup of order 1 of Elliptic Curve defined by y^2 = x^3 + 2815805388*x + 94443526967868 over Rational Field
C1.torsion_subgroup()
Torsion subgroup of order 1 of Elliptic Curve defined by y^2 = x^3 + 2265722465761*x + 3154189403034549278 over Rational Field
C5.torsion_subgroup()
Torsion subgroup of order 1 of Elliptic Curve defined by y^2 = x^3 + 2639323244332897*x + 20156152630838819347102 over Rational Field
C7.torsion_subgroup()
Torsion subgroup of order 1 of Elliptic Curve defined by y^2 = x^3 + 349942184229228*x + 11167797929528591502588 over Rational Field
C9.torsion_subgroup()
Torsion subgroup of order 1 of Elliptic Curve defined by y^2 = x^3 + 9243195710310751148*x + 761969307339454319105751548 over Rational Field
"""

"""
Rational ponts
E.rational_points(bound=1000): # 'bound' is naive height
What do they mean?
"""