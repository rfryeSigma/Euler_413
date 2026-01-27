"""
Database of known solutions
"""
from gmpy2 import digits, is_prime
from math import ceil, prod

# dictionary of 20 known and verified solutions.
# key is {1st letter of info}, {5 digits of d}, '_', 
#   {number of digits in d}
# abcd is sorted tuple (a, b, c, d) in a^4 + b^4 + c^4 = d^4
# tuvw tuple (t=a//8, u=b//8, v=(d-c)//2, w=(d+c)//2) where
#    a and b are chosen as even terms, not smallest terms.
# v_evens and v_primes are the factors of v.
# w_evens and w_primes are the factors of w.

known = dict(
F42248_6=dict(info='Frye (1988)',
	abcd=(95_800, 217_519, 414_560, 422_481),
	tuvw=(11_975, 51_820, 102_481, 320_000),
	v_evens=2**0, v_primes={102_481: 1},
	w_evens=2**9, w_primes={5: 4},
	common={5: 1},
	factors={102_481: 1, 26_777: 1, 4_216_393: 1},
	others={},
),
M28130_7=dict(info='MacLeod (1997)',
	abcd=(673_865, 1_390_400, 2_767_624, 2_813_001),
	tuvw=(173_800, 345_953, 1_069_568, 1_743_433),
	v_evens=2**9, v_primes={2_089: 1},
	w_evens=2**0, w_primes={1_743_433: 1},
	common={},
	factors={73: 1, 2_089: 1, 1_743_433: 1, 57_308_689_481: 1},
	others={},
),
B87074_7=dict(info='Bernstein (2001)',
	abcd=(1_705_575, 5_507_880, 8_332_208, 8_707_481),
	tuvw=(688_485, 1_041_526, 3_500_953, 5_206_528),
	v_evens=2**0, v_primes={3_500_953: 1},
	w_evens=2**9, w_primes={10_169: 1},
	common={},
	factors={13_031_561: 1, 17: 1, 1_297: 1, 137: 1, 10_169: 1, 3_500_953: 1},
	others={},
),
T12197_8=dict(info='Tomita (Curve 4)',
	abcd=(5_870_000, 8_282_543, 11_289_040, 12_197_457),
	tuvw=(733_750, 1_411_130, 1_957_457, 10_240_000),
	v_evens=2**0, v_primes={601: 1, 3_257: 1},
	w_evens=2**14, w_primes={5: 4},
	common={2: 1, 5: 1},
	factors={2: 1, 19_213_229_257: 1, 601: 1, 3_257: 1, 5_657: 1},
	others={},
),
B16003_8=dict(info='Bernstein (2001)',
	abcd=(4_479_031, 12_552_200, 14_173_720, 16_003_017),
	tuvw=(1_569_025, 1_771_715, 5_761_993, 10_241_024),
	v_evens=2**0, v_primes={5_761_993: 1},
	w_evens=2**10, w_primes={73: 1, 137: 1},
	common={5: 1},
	factors={2: 1, 5_761_993: 1, 220_926_617_441: 1, 73: 1, 137: 1},
	others={},
),
B16430_8=dict(info='Bernstein (2001)',
	abcd=(3_642_840, 7_028_600, 16_281_009, 16_430_513),
	tuvw=(455_355, 878_575, 74_752, 16_355_761),
	v_evens=2**10, v_primes={73: 1},
	w_evens=2**0, w_primes={41: 1, 398_921: 1},
	common={5: 1},
	factors={2: 1, 398_921: 1, 30_039_049: 1, 41: 1, 73: 1, 14_249: 1},
	others={},
),
E20615_8=dict(info='Elkies (1998)',
	abcd=(2_682_440, 15_365_639, 18_796_760, 20_615_673),
	tuvw=(335_305, 2_349_595, 2_625_017, 17_990_656),
	v_evens=2**0, v_primes={2_625_017: 1},
	w_evens=2**10, w_primes={17_569: 1},
	common={5: 1},
	factors={2: 1, 2_625_017: 1, 97: 1, 353: 1, 761: 1, 17_569: 1, 20_297: 1},
	others={},
),
G44310_8=dict(info='Gerbicz (2006)',
	abcd=(2_164_632, 31_669_120, 41_084_175, 44_310_257),
	tuvw=(270_579, 3_958_640, 1_613_041, 42_697_216),
	v_evens=2**0, v_primes={1_613_041: 1},
	w_evens=2**9, w_primes={89: 1, 937: 1},
	common={},
	factors={9_004_361: 1, 11_926_601: 1, 17: 1, 89: 1, 937: 1, 1_613_041: 1},
	others={},
),
G68711_8=dict(info='Gerbicz (2006)',
	abcd=(10_409_096, 42_878_560, 65_932_985, 68_711_097),
	tuvw=(1_301_137, 5_359_820, 1_389_056, 67_322_041),
	v_evens=2**9, v_primes={2_713: 1},
	w_evens=2**0, w_primes={41: 1, 457: 1, 3_593: 1},
	common={},
	factors={3_593: 1, 626_177: 1, 7_241_062_321: 1, 41: 1, 457: 1, 2_713: 1},
	others={},
),
G11711_9=dict(info='Gerbicz (2006)',
	abcd=(34_918_520, 87_865_617, 106_161_120, 117_112_081),
	tuvw=(4_364_815, 13_270_140, 14_623_232, 102_488_849),
	v_evens=2**9, v_primes={13: 4},
	w_evens=2**0, w_primes={281: 1, 569: 1, 641: 1},
	common={5: 1, 13: 1, 53: 1},
	factors={641: 1, 569: 1, 6_449: 1, 281: 1, 337: 1},
	others={},
),
R14508_9=dict(info='Rathmann (2006)',
	abcd=(1_841_160, 121_952_168, 122_055_375, 145_087_793),
	tuvw=(230_145, 15_244_021, 11_516_209, 133_571_584),
	v_evens=2**0, v_primes={313: 1, 36_793: 1},
	w_evens=2**10, w_primes={17: 1, 7_673: 1},
	common={},
	factors={2: 1, 2_752_947_024_353: 1, 17: 1, 313: 1, 6_529: 1, 7_673: 1, 36_793: 1},
	others={},
),
R15664_9=dict(info='Rathmann (2006)',
	abcd=(27_450_160, 108_644_015, 146_627_384, 156_646_737),
	tuvw=(3_431_270, 18_328_423, 24_001_361, 132_645_376),
	v_evens=2**0, v_primes={24_001_361: 1},
	w_evens=2**9, w_primes={449: 1, 577: 1},
	common={},
	factors={62_576_009: 1, 161_233: 1, 24_001_361: 1, 449: 1, 577: 1, 1_801: 1},
	others={},
),
T58984_9=dict(info='Tomita (2006)',
	abcd=(186_668_000, 260_052_385, 582_665_296, 589_845_921),
	tuvw=(23_333_500, 72_833_162, 164_896_768, 424_949_153),
	v_evens=2**13, v_primes={20_129: 1},
	w_evens=2**0, w_primes={17: 1, 24_997_009: 1},
	common={2: 1},
	factors={24_997_009: 1, 17: 1, 41: 2, 233: 1, 20_129: 1, 54_570_001: 1, 9_721: 1},
	others={},
),
M63852_9=dict(info='MacLeod (1997)',
	abcd=(219_076_465, 275_156_240, 630_662_624, 638_523_249),
	tuvw=(34_394_530, 78_832_828, 209_723_392, 428_799_857),
	v_evens=2**13, v_primes={25_601: 1},
	w_evens=2**0, w_primes={17: 1, 113: 1, 223_217: 1},
	common={2: 1},
	factors={223_217: 1, 17: 1, 113: 1, 1_697: 1, 25_601: 1, 3_341_489: 1, 40_182_161: 1},
	others={},
),
G87382_9=dict(info='Gerbicz (2009)',
	abcd=(558_424_440, 606_710_871, 769_321_280, 873_822_121),
	tuvw=(69_803_055, 96_165_160, 133_555_625, 740_266_496),
	v_evens=2**0, v_primes={5: 4, 7: 4, 89: 1},
	w_evens=2**9, w_primes={17: 1, 85_049: 1},
	common={5: 1, 7: 1},
	factors={22_769: 1, 17: 1, 89: 1, 16_273: 1, 85_049: 1, 1_527_128_593: 1},
	others={},
),
G12597_10=dict(info='Gerbicz (2009)',
	abcd=(588_903_336, 859_396_455, 1_166_705_840, 1_259_768_473),
	tuvw=(73_612_917, 145_838_230, 200_186_009, 1_059_582_464),
	v_evens=2**0, v_primes={89: 1, 2_249_281: 1},
	w_evens=2**9, w_primes={2_069_497: 1},
	common={41: 1},
	factors={431_790_169: 1, 89: 1, 953: 1, 2_069_497: 1, 2_249_281: 1},
	others={},
),
T16791_10=dict(info='Tomita (2006)',
	abcd=(50_237_800, 632_671_960, 1_670_617_271, 1_679_142_729),
	tuvw=(6_279_725, 79_083_995, 4_262_729, 1_674_880_000),
	v_evens=2**0, v_primes={41: 1, 103_969: 1},
	w_evens=2**10, w_primes={5: 4, 2_617: 1},
	common={5: 1},
	factors={2: 1, 17: 1, 41: 1, 103_969: 1, 2_617: 1, 81_049: 1, 34_361: 1, 59_252_657: 1},
	others={},
),
G17878_10=dict(info='Gerbicz (2009)',
	abcd=(686_398_000, 1_237_796_960, 1_662_997_663, 1_787_882_337),
	tuvw=(85_799_750, 154_724_620, 62_442_337, 1_725_440_000),
	v_evens=2**0, v_primes={1_097: 1, 56_921: 1},
	w_evens=2**13, w_primes={5: 4, 337: 1},
	common={2: 1, 5: 1},
	factors={5_297: 1, 56_921: 1, 84_526_369: 1, 17: 1, 337: 1, 457: 1, 857: 1, 1_097: 1},
	others={},
),
G18717_10=dict(info='Gerbicz (2009)',
	abcd=(92_622_401, 1_553_556_440, 1_593_513_080, 1_871_713_857),
	tuvw=(194_194_555, 199_189_135, 889_545_728, 982_168_129),
	v_evens=2**10, v_primes={868_697: 1},
	w_evens=2**0, w_primes={982_168_129: 1},
	common={5: 1},
	factors={2: 1, 17: 1, 3_729_457: 1, 868_697: 1},
	others={43_523_359_439_352_337: 1},
),
T29999_14=dict(info='Tomita (2006)',
	abcd=(7_592_431_981_391, 22_495_595_284_040, 27_239_791_692_640, 29_999_857_938_609),
	tuvw=(2_811_949_410_505, 3_404_973_961_580, 11_203_712_978_609, 18_796_144_960_000),
	v_evens=2**0, v_primes={11_203_712_978_609: 1},
	w_evens=2**9, w_primes={5: 4, 41: 1, 89: 1, 16_097: 1},
	common={5: 1},
	factors={16_097: 1, 41: 1, 73: 1, 89: 1, 12_329: 1},
	others={5_960_489_902_302_355_385_651_356_683_024_737: 1},
),
)

def make_known(verbose: bool=True):
    """
    Verify known table.
    Report or add a potential new field.
    If v: copy new table entries.
    """
    for key, val in known.items():
        # Verify entry
        info = val['info']
        a, b, c, d = abcd = val['abcd']
        verify_key(key, info, d)
        verify_abcd(abcd, False)

        t, u, v, w = tuvw = val['tuvw']
        verify_tuvw(abcd, tuvw)

        v_evens = val['v_evens']
        v_primes = val['v_primes']
        for k in v_primes.keys():
            assert is_prime(k)
        assert v == v_evens * prod([k**e for k, e in v_primes.items()])

        w_evens = val['w_evens']
        w_primes = val['w_primes']
        for k in w_primes.keys():
            assert is_prime(k)
        assert w == w_evens * prod([k**e for k, e in w_primes.items()])

        common = val['common']
        factors = val['factors']
        others = val['others']

        # Code to add to an entry
        pass

        # Opional output modified entry
        if not verbose: continue
        vz = (v & -v).bit_length() - 1
        vp = (f'{k:_}: {e}' for k, e in v_primes.items())
        wz = (w & -w).bit_length() - 1
        wp =[f'{k:_}: {e}' for k, e in w_primes.items()]
        cp = (f'{int(k):_}: {e}' for k, e in common.items())
        fp = (f'{int(k):_}: {e}' for k, e in factors.items())
        of = (f'{int(k):_}: {e}' for k, e in others.items())
        print(f"{key}=dict(info='{info}',"
              f"\n\tabcd=({a:_}, {b:_}, {c:_}, {d:_}),"
              f"\n\ttuvw=({t:_}, {u:_}, {v:_}, {w:_}),"
              f"\n\tv_evens=2**{vz}, v_primes={{{', '.join(map(str, vp))}}},"
              f"\n\tw_evens=2**{wz}, w_primes={{{', '.join(map(str, wp))}}},"
              f"\n\tcommon={{{', '.join(map(str, cp))}}},"
              f"\n\tfactors={{{', '.join(map(str, fp))}}},"
              f"\n\tothers={{{', '.join(map(str, of))}}},"
              "\n),")

# dictionary of known, but too long to use solutions
too_long = { 
    'E91617.71': dict(info='Elkies 1988',
                     abcd=
          (1439965710648954492268506771833175267850201426615300442218292336336633, 
           4417264698994538496943597489754952845854672497179047898864124209346920, 
           9033964577482532388059482429398457291004947925005743028147465732645880, 
           9161781830035436847832452398267266038227002962257243662070370888722169)),
}
"""
1.37 10^15,Large Quartet,Found by Elkies (Rank 3 curve)
2.02 10^18,Large Quartet,Bremner
2.77 10^18,Large Integer,Large Integer,Large Integer,Bremner (2015)
"""

# notes on supposed solutions that fail tests
failures = dict( 
) # discarded many supposed solutions that didn't get past modulo 1000.

def verify_tuvw(abcd, tuvw) -> bool:
    """Check tuvw matches abcd
    """
    a, b, c, d = abcd
    t, u, v, w = tuvw
    assert t < u and v < w
    assert v % 2**9 == 0 or w % 2**9 == 0
    cc = next(x for x in abcd[:-1] if x % 2 > 0)
    assert d == w + v and cc == w - v
    M = d**4 - cc**4
    assert t**4 + u**4 == M // 2**12
    assert M // 2**3 == v * w**3 + w * v**3
    return True

def verify_abcd(solution: tuple, verbose: bool=True) -> bool:
    """Check solution a^4 + b^4 + c^4 == d^4.
    """
    a, b, c, d = solution
    assert a <= b <= c <= d
    lhs = a ** 4 + b ** 4 + c ** 4
    rhs = d ** 4
    if lhs == rhs:
        if verbose: print(f'Match {int(rhs):_}')
        return True
    else:
        if verbose: print(f'FAIL\n{int(lhs):_}\n{int(rhs):_}')
        return False

def verify_abcd_modulo_10(solution: tuple) -> int:
    """Check solution at successive powers of 10.
    """
    a, b, c, d = solution
    assert a <= b <= c <= d
    d4 = d ** 4
    e = 1
    modulus = 10
    while True:
        lhs = (pow(a, 4, modulus) + pow(b, 4, modulus) + pow(c, 4, modulus)) % modulus
        rhs = pow(d, 4, modulus)
        print(e, modulus, lhs, rhs)
        if lhs != rhs:
            print('fail')
            return e
        if rhs >= d4:
          print('All powers pass')
          return e
        e += 1
        modulus *= 10

def generate_key(info: str, d: int) -> str:
    """Generate dict key for a solution.
    The key descibes the entry info and d as "Lxyz_n, 
    where L is the first letter of info,
    xyz are first five digits of d,
    and n is number of digits of d"
    """
    L = info[0]
    assert L.isalpha()
    digs = digits(d)
    key = L + digs[:5] + '_' + str(len(digs))
    return key

def verify_key(key: str, info: str, d: int) -> bool:
    if not info.startswith(key[0]):
        import pdb; pdb.set_trace()
    assert info.startswith(key[0])
    digs, num = str(key[1:]).split('_')
    n = eval(num)
    assert eval(digs) == d // 10**(n - 5)
    return True

# This is how p_max is estimated:
# p_max is an estimate of the max prime needed for solve.
# t^4 + b^4 = m; 2^9 * m = w * v^3  + v * w^3 = x * (y^3 + x^2 * y)
# p_max is when x is w and is prime.
# [ceil(w / v) for val in solutions.known.values() for v, w in [val['tuvw'][2:]]]
#  = [4, 2, 2, 6, 2, 219, 7, 27, 49, 8, 12, 6, 3, 3, 6, 6, 393, 28, 2, 2]
# Select r = 2 as an initial low estimate estimate of the ratio.
# 2^9 * m = w * v^3 + v w^3 ~= w^4 * (1/r^3 + 1/r)
# The solver builds LHS of x starting with 2^9.
# 2^9 * m > (2^9 * p_max)^4 * (1/r^3 + 1/r)
# p_max < 2*-6 * (m * (r^3)/(1 + r^2) * 2^-3)^(1/4)
def estimate_max_prime():
    """Estimate largest prime needed to factor known solutions.
    """
    for key, val in known.items():
        # Extract entry
        t, u, v, w = tuvw = val['tuvw']
        v_evens = val['v_evens']
        v_primes = val['v_primes']
        w_evens = val['w_evens']
        w_primes = val['w_primes']

        # Calculate p_max
        m = t**4 + u**4
        r = 2
        # p_max < 2*-6 * (m * (r^3)/(1 + r^2) * 2^-3)^(1/4)
        p_max = ceil(2**-6 * (m * r**3 * (1 + r*r)**-1)**(1/4))
        print(key, f'{p_max:_}'
              f'\n\tv: {v_evens}, {v_primes}'
              f'\n\tw: {w_evens}, {w_primes}')
"""
This is the output for the 10 and 14 digit solutions with comments
G12597_10 2_603_455 # Just enough for the worst case
	v: 1, {89: 1, 2249281: 1}
	w: 512, {2069497: 1}
T16791_10 1_389_771
	v: 1, {41: 1, 103969: 1}
	w: 1024, {5: 4, 2617: 1}
G17878_10 2_781_118
	v: 1, {1097: 1, 56921: 1}
	w: 8192, {5: 4, 337: 1}
G18717_10 4_111_487
	v: 1024, {868697: 1}
	w: 1, {982168129: 1}
T29999_14 65_831_370_305
	v: 1, {11203712978609: 1}
	w: 512, {5: 4, 41: 1, 89: 1, 16097: 1}

tree1_18 has primes 17 thru 7_766_713
Maximum prime in 20 indexes at level 12 is
>>> factoring.tree1[0][20 * 2**12 - 1]
4_686_281
which is more than the penultimate estimate.

The final T29999_14 estimate is not practical,
but the factors needed all happen to be small.
"""

def solve_all():
    """Solve all known solutions.
    """
    from logging_413 import V, IntFlag
    from factoring import factor_tu
    from solving import solve_factors
    from time import time
    start = time()
    for key, val in known.items():
        t, u, v, w = tuvw = val['tuvw']
        print(f'{key}: tuvw=({t:_}, {u:_}, {v:_}, {w:_}),')
 
        common, factors, others = factor_tu(t, u, vV=V.NONE)
        cp = (f'{int(k):_}: {e}' for k, e in common.items())
        fp = (f'{int(k):_}: {e}' for k, e in factors.items())
        of = (f'{int(k):_}: {e}' for k, e in others.items())
        print(f"\tcommon={{{', '.join(map(str, cp))}}},")
        print(f"\tfactors={{{', '.join(map(str, fp))}}},")
        print(f"\tothers={{{', '.join(map(str, of))}}},")

        sol = solve_factors(common, factors, others, vV=V.NONE)
        if sol is None:
            print('\tNo solution')
        else:
            s1, s2 = map(int, sol)
            print(f'\tsolution: ({s1:_}, {s2:_})')
    elapsed = time() - start
    print(f'Elapsed: {elapsed:.6f}s')
"""
This is an extract of the output:
F42248_6: tuvw=(11_975, 51_820, 102_481, 320_000),
	common={5: 1},
	factors={102_481: 1, 26_777: 1, 4_216_393: 1},
	others={},
	solution: (102_481, 320_000)
...
    T29999_14: tuvw=(2_811_949_410_505, 3_404_973_961_580, 11_203_712_978_609, 18_796_144_960_000),
	common={5: 1},
	factors={16_097: 1, 41: 1, 73: 1, 89: 1, 12_329: 1},
	others={5_960_489_902_302_355_385_651_356_683_024_737: 1},
	solution: (11_203_712_978_609, 18_796_144_960_000)
Elapsed: 0.008516s
"""
