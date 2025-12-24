# Decompose known solutions to a^4 + b^4 + c^4 = d^4.
import math
import sys
from factoring_strategy import factoring_strategy

known_solutions = { # key format: (a,b,c,d)
    (95800, 217519, 414560, 422481): dict(
        index=0, # solution #0: Roger Frye 1988
        tuvw=(11975, 51820, 102481, 320000),
        t_evens=1, t_primes={5: 2, 479: 1},
        u_evens=4, u_primes={5: 1, 2591: 1},
        v_evens=1, v_primes={102481: 1},
        w_evens=512, w_primes={5: 4},
        v2_w2=112902355361, v2_w2_evens=1, v2_w2_primes={26777: 1, 4216393: 1},
    ),

    (2682440, 15365639, 18796760, 20615673): dict( # This solution has the 2nd smallest a.
        index=1, # 1: Noam Elkies 1988
        tuvw=(335305, 2349595, 2625017, 17990656),
        t_evens=1, t_primes={5: 1, 67061: 1},
        u_evens=1, u_primes={5: 1, 469919: 1},
        v_evens=1, v_primes={2625017: 1},
        w_evens=1024, w_primes={17569: 1},
        v2_w2=330554417560625, v2_w2_evens=1, v2_w2_primes={5: 4, 97: 1, 353: 1, 761: 1, 20297: 1},
    ),

    (34918520, 87865617, 106161120, 117112081): dict(
        index=2, # 2: Robert Gerbicz 2006
        tuvw=(4364815, 13270140, 14623232, 102488849),
        t_evens=1, t_primes={5: 1, 7: 1, 13: 1, 53: 1, 181: 1},
        u_evens=4, u_primes={3: 2, 5: 1, 13: 1, 53: 1, 107: 1},
        v_evens=512, v_primes={13: 4},
        w_evens=1, w_primes={281: 1, 569: 1, 641: 1},
        v2_w2=10717803083470625, v2_w2_evens=1, v2_w2_primes={5: 4, 53: 4, 337: 1, 6449: 1},
    ),

    # This solution shows that there may be no a,b,c without factors of 2 or 5
    # Also d^4-c^4 has 2^16 factor, 4 more than usual
    (186668000, 260052385, 582665296, 589845921): dict(
        index=3, # 3: Seiji Tomita 2006
        tuvw=(23333500, 72833162, 164896768, 424949153),
        t_evens=4, t_primes={5: 3, 23: 1, 2029: 1},
        u_evens=2, u_primes={36416581: 1},
        v_evens=8192, v_primes={20129: 1},
        w_evens=1, w_primes={17: 1, 24997009: 1},
        v2_w2=207772726732263233 , v2_w2_evens=1, v2_w2_primes={41: 2, 233: 1, 9721: 1, 54570001: 1},
     ),

    (50237800, 632671960, 1670617271, 1679142729): dict(
        index=4, # 4: Seiji Tomita 2006
        tuvw=(6279725, 79083995, 4262729, 1674880000),
        t_evens=1, t_primes={5: 2, 239: 1, 1051: 1},
        u_evens=1, u_primes={5: 1, 15816799: 1},
        v_evens=1, v_primes={41: 1, 103969: 1},
        w_evens=1024, w_primes={5: 4, 2617: 1},
        v2_w2=2805241185258527441,v2_w2_evens=1, v2_w2_primes={17: 1, 34361: 1, 81049: 1, 59252657: 1},
    ),

    (7592431981391, 22495595284040, 27239791692640, 29999857938609): dict(
        index=5, # 5: Seiji Tomita 2006
        tuvw=(2811949410505, 3404973961580, 11203712978609, 18796144960000),
        t_evens=1, t_primes={5: 1, 19: 1, 29599467479: 1},
        u_evens=4, u_primes={5: 1, 170248698079: 1},
        v_evens=1, v_primes={11203712978609: 1},
        w_evens=512, w_primes={5: 4, 41: 1, 89: 1, 16097: 1},
        v2_w2=478818249864385152491574881, v2_w2_evens=1,
        v2_w2_primes={73: 1, 12329: 1, 532010228544999874993: 1},
    ),

    (558424440, 606710871, 769321280, 873822121): dict(
        index=6, # 6: Robert Gerbicz 2009
        tuvw=(69803055, 96165160, 133555625, 740266496),
        t_evens=1, t_prime={3: 2, 5: 1, 7: 1, 19: 1, 107: 1, 109: 1},
        u_evens=8, u_primes={5: 1, 7: 1, 13: 1, 29: 1, 911: 1},
        v_evens=1, v_primes={5: 4, 7: 4, 89: 1},
        w_evens=512, w_primes={17: 1, 85049: 1},
        v2_w2=565831590069258641, v2_w2_evens=1, v2_w2_primes={16273: 1, 22769: 1, 1527128593: 1},
    ),

    (588903336, 859396455, 1166705840, 1259768473): dict(
        index=7, # 7: Robert Gerbicz 2009
        tuvw=(73612917, 145838230, 200186009, 1059582464),
        t_evens=1, t_primes={3: 2, 7: 1, 41: 1, 28499: 1},
        u_evens=2, u_primes={5: 1, 41: 1, 67: 1, 5309: 1},
        v_evens=1, v_primes={89: 1, 2249281: 1},
        w_evens=512, w_primes={2069497: 1},
        v2_w2=1162789436215659377, v2_w2_evens=1, v2_w2_primes={41: 4, 953: 1, 431790169: 1},
    ),

    # d^4-c^4 has 2^16 factor, 4 more than usual
    (686398000, 1237796960, 1662997663, 1787882337): dict(
        index=8, # 8: Robert Gerbicz 2009
        tuvw=(85799750, 154724620, 62442337, 1725440000),
        t_evens=2, t_primes={5: 3, 343199: 1}, 
        u_evens=4, u_primes={5: 1, 71: 1, 108961: 1},
        v_evens=1, v_primes={1097: 1, 56921: 1},
        w_evens=8192, w_primes={5: 4, 337: 1},
        v2_w2=2981042239050021569, v2_w2_evens=1,
	    v2_w2_primes={17: 1, 457: 1, 857: 1, 5297: 1, 84526369: 1},
	    ),

    (92622401, 1553556440, 1593513080, 1871713857): dict(
        index=9, # 9: Robert Gerbicz 2009
        tuvw=(194194555, 199189135, 889545728, 982168129),
        t_evens=1, t_primes={5: 1, 38838911: 1},
        u_evens=1, u_primes={5: 1, 257: 1, 379: 1, 409: 1},
        v_evens=1024, v_primes={868697: 1},
        w_evens=1, w_primes={982168129: 1},
        v2_w2=1755945835826410625, v2_w2_evens=1,
        v2_w2_primes={5: 4, 17: 1, 3729457: 1, 44313553: 1},
    ),

    (1439965710648954492268506771833175267850201426615300442218292336336633, 
     4417264698994538496943597489754952845854672497179047898864124209346920, 
     9033964577482532388059482429398457291004947925005743028147465732645880, 
     9161781830035436847832452398267266038227002962257243662070370888722169): dict(
        index=10, # Also 10: Noam Elkies 1988, 71 digit solution
        tuvw=(552158087374317312117949686219369105731834062147380987358015526168365, 
              1129245572185316548507435303674807161375618490625717878518433216580735, 
              3860908059693241177781972813217045385188400767820971609926039276192768, 
              5300873770342195670050479585050220653038602194436272052144331612529401),
    # a = 1439965710648954492268506771833175267850201426615300442218292336336633 
    #   factors: {83: 1, 26711: 1, 15365639: 1, 42270100962827684990671492004308319650144870955038472419: 1}
    # b = 4417264698994538496943597489754952845854672497179047898864124209346920
    #   factors {2: 3, 5:1, 542603: 1} [203521944174402763021195860037400864253177391996498724613765691,]
    # c = 9033964577482532388059482429398457291004947925005743028147465732645880
    #   factors {2: 3, 5: 1, 739: 1, 13183: 1, 295879: 1, 78351183214844218506223148101631054732991620013524560889: 1
    # d = 9161781830035436847832452398267266038227002962257243662070370888722169
    #   factors {3: 1, 13: 2, 173: 1, 3433: 1, 6007: 1, 5065174920412934124179876558822872718782837048261723640809: 1}
    #
    # x = a = 1439965710648954492268506771833175267850201426615300442218292336336633:
    # d-x = 7721816119386482355563945626434090770376801535641943219852078552385536 
    #   factors: {2: 11, 2633: 1, 4409: 1, 324786930837975462929657126375337026762581115572406886530081: 1}
    # d+x = 10601747540684391340100959170100441306077204388872544104288663225058802 
    #   factors: {2: 1, 137: 1, 38692509272570771314237077263140296737508045214863299650688551916273: 1}
    # d^2+x^2 = 86011747549012226770024291965474831199593742806467880998972755694111066404078098578470531153488908905607830781751532005752181596156369841250 
    #   factors: {2: 1, 5: 4, 17: 1, 929: 1, 38329: 1, 527489: 1, 215497495402636451793613401133356993585929216173254203034205207089732330345050773904457735361825500277702818296492875280281: 1}
    # I later decomposed these into primes with help from WolframAlpha.
     ),
}

def verify_keys():
    """Verify that the keys to known_solutions satisfy a^4 + b^4 + c^4 = d^4.
    """
    verified = True
    for sol, info in known_solutions.items():
        a, b, c, d = sol
        assert a < b < c < d, 'Expected a < b < c < d in solution {}'.format(sol)
        lhs = a**4 + b**4 + c**4
        rhs = d**4
        if lhs != rhs:
            verified = False
            print(f'Solution #{info['index']} {sol} FAILED verification: {lhs} != {rhs}')
        else:
            print(f'Solution #{info['index']} {sol} verified')
    assert verified, 'Some solutions failed verification.'

def redefine_abcd():
    """Redefine a, b, c, d solution in terms of t, u, v, w
    """
    import pdb; pdb.set_trace()
    verified = True
    for sol, info in known_solutions.items():
        a, b, c, d = sol
        assert d % 2 == 1
        if a % 8 == 0:
            t = a // 8
            if b % 8 == 0:
                u = b // 8
                assert c % 2 == 1
                v = (d - c) // 2
                w = (d + c) // 2
                assert v * w * 4 == (d - c) * (d + c)
                assert v + w == d
                assert w - v == c
                M = a**4 + b**4
                M2 = d**4 - c**4
                assert M == M2
                assert 2**12 * (t**4 + u**4) == M
                assert 2**3 * w * v * (w**2 + v**2) == M
            else:
                assert c % 8 == 0
                u = c // 8
                assert b % 2 == 1
                v = (d - b) // 2
                w = (d + b) // 2
                assert v * w * 4 == (d - b) * (d + b)
                assert v + w == d
                assert w - v == b
                M = a**4 + c**4
                M2 = d**4 - b**4
                assert M == M2
                assert 2**12 * (t**4 + u**4) == M
                assert 2**3 * w * v * (w**2 + v**2) == M
        else:
            assert a % 2 == 1 and b % 8 == 0 and c % 8 == 0
            t = b // 8
            u = c // 8
            v = (d - a) // 2
            w = (d + a) // 2
            assert v * w * 4 == (d - a) * (d + a)
            assert v + w == d
            assert w - v == a
            M = b**4 + c**4
            M2 = d**4 - a**4
            assert M == M2
            assert 2**12 * (t**4 + u**4) == M
            assert 2**3 * w * v * (w**2 + v**2) == M
        tuvw = (t, u, v, w)
        lhs = t**4 + u**4
        rhs = v * w * (v**2 + w**2)
        if lhs * 2**9 != rhs:
            verified = False
            print(f'Solution #{info['index']} {tuvw} FAILED: {rhs} / {lhs} = {rhs / lhs}')
        else:
            print(f'Solution #{info['index']} {tuvw} verified')
    assert verified, 'Some solutions failed verification.'

def verify_tuvw():
    """Verify that the tuvw entry in known_solutions satisfy equation.
    """
    verified = True
    for info in known_solutions.values():
        t, u, v, w = tuvw = info['tuvw']
        lhs = t**4 + u**4
        rhs = v * w * (v**2 + w**2)
        if lhs * 2**9 != rhs:
            verified = False
            print(f'Solution #{info['index']} {tuvw} FAILED: {rhs} / {lhs} = {rhs / lhs}')
        else:
            print(f'Solution #{info['index']} {tuvw} verified')
    assert verified, 'Some solutions failed verification.'

def factor_tuvw(inx: int, verbose: bool=False,):
    """Factor the tuvd parameters for one known solution
    Args:
        inx (int): Index of known solution to factor.
    """
    assert 0 <= inx < len(known_solutions), 'Invalid solution index {}'.format(inx)
    for info in known_solutions.values():
        if info['index'] == inx:
            break
    else:
        print(f'Failed to find solution with index {inx}')
        return 1
    t, u, v, w = tuvw = info['tuvw'] 
    print(f'Factoring tuvw {tuvw}...')

    evens, factors, others = factoring_strategy(t, verbose=verbose)
    print(f't: {t} -> evens: {evens},'
              f'\n\tprimes: {factors},\n\tothers: {others}')

    evens, factors, others = factoring_strategy(u, verbose=verbose)
    print(f'u: {u} -> evens: {evens},'
              f'\n\tprimes: {factors},\n\tothers: {others}')

    evens, factors, others = factoring_strategy(v, verbose=verbose)
    print(f'v: {v} -> evens: {evens},'
              f'\n\tprimes: {factors},\n\tothers: {others}')

    evens, factors, others = factoring_strategy(w, verbose=verbose)
    print(f'w: {w} -> evens: {evens},'
              f'\n\tprimes: {factors},\n\tothers: {others}')

    v2w2 = v**2 + w**2
    evens, factors, others = factoring_strategy(v2w2, verbose=verbose)
    print(f'v^2+w^2: {v2w2} -> evens: {evens},'
              f'\n\tprimes: {factors},\n\tothers: {others}')

def pmod8():
    """Verify that prime factors of v, w, v^2+w^2 are all 1 mod 8 unless exponent 4.
    """
    verified = True
    pmod8_vals = {1:0, 3:0, 5:0, 7:0}
    for info in known_solutions.values():
        index = info['index']
        if index == 10: # Except #10 Elkies 71 digit
            continue
        for key in ('v_primes', 'w_primes', 'v2_w2_primes'):
            for p, e in info[key].items():
                val = p % 8
                pmod8_vals[val] += e
                if val != 1 and e != 4:
                    verified = False
                    print(f'prime {p} % 8 = {val} in {key} in solution {index}')
    print(f'pmod8_vals {pmod8_vals}')
    assert verified, 'Some solutions failed verification.'

def residues():
    """Report a dict of all v, w, v^2+w^2 residues mod 5, 17, 29
    This would be to verify LPS rules for vw_search.
    """
    print('Collected residues is not yet implemented.')

def fit_parameterizations():
    """Fit parameterizations to known solutions.
    """
    # Can I identify known solutions in Elkies' parameterization?
    print('Fitting parameterizations is not yet implemented.')

## Main
# ----------------------------------------------------------------------------
import argparse
import sys

def main(argv=None):
    """Command-line dispatcher for `decompose_known_solutions.py`.
    """
    parser = argparse.ArgumentParser(prog='decompose_known_solutions',
        description='Decompose known solutions to a^4 + b^4 + c^4 = d^4.')
    subparsers = parser.add_subparsers(dest='command')  

    # verify_keys command
    p_verify_keys = subparsers.add_parser('verify_keys', 
        help='verify that the keys are solutions')

    # redefine_abcd command
    p_edefine_abcd = subparsers.add_parser('redefine_abcd', 
        help='Redefine a, b, c, d solution in terms of t, u, v, w')

    # verify_tuvw command
    p_verify_tuvw = subparsers.add_parser('verify_tuvw', 
        help='Verify that the tuvw entrys satisfy equation')

    # factor command
    p_factor_tuvw = subparsers.add_parser('factor_tuvw', 
        help='Factor the tuvd parameters for one known solution')
    p_factor_tuvw.add_argument('index', type=int, 
        help='index to identify which solution to factor')
    p_factor_tuvw.add_argument('-v', action='store_true', 
        help='whether to log progress')

    # pmod8 command
    p_pmod8 = subparsers.add_parser('pmod8', 
        help='verify that prime factors are 1 mod 8')
   
    args, rest = parser.parse_known_args(argv)

    if args.command == 'verify_keys':
        verify_keys()
        return 0

    if args.command == 'redefine_abcd':
        redefine_abcd()
        return 0

    if args.command == 'verify_tuvw':
        verify_tuvw()
        return 0

    if args.command == 'factor_tuvw':
        factor_tuvw(args.index, verbose=args.v)
        return 0

    if args.command == 'pmod8':
        pmod8()
        return 0

    parser.error('unknown command')

if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
