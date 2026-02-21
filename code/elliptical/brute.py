"""
Run a brute force search over small rationals for solutions using Pieza's pairs method.
https://math.stackexchange.com/questions/4857229/on-why-solutions-to-x4y4z4-1-come-in-pairs

"""
from solutions import known
import csv
from datetime import datetime
from itertools import combinations
from math import isqrt, gcd
from multiprocessing import Process, Queue
from numpy import searchsorted
from sage.all import Integer, QQ, Rational, lcm, numerator, denominator
from timeit import repeat, time

def abcd_to_xyz(abcd: tuple) -> tuple:
    """Convert (a, b, c, d) representing a^4 + b^4 + c^4 = d^4 
    to rational form (x, y, z) representing x^4 + y^4 + z^4 = 1
    """
    A, B, C, D = map(QQ, abcd)
    x, y, z = A/D, B/D, C/D
    return x, y, z

def xyz_to_u(x :Rational, y :Rational, z :Rational) -> Rational:
    """Convert (x, y, z) representing x^4 + y^4 + z^4 = 1
    to single parameter u.
    """
    x_y = x - y
    u = (x_y*x_y - z*z - 1) / (x*x - x*y + y*y + x_y)
    return u

def abcd_to_u_set(abcd: tuple) -> set:
    """Convert (a, b, c, d) to (x, y, z) and then consider
    3 permutations and 8 combinations of signs on (x, y, z).
    Since the u, v, w formulas rotate (x, y, z), only swap y, z.
    to generate the set of all EC parameters u.
    There will be 12 duplicates, so return unique 12.
    """
    u_set = set()
    xx, yy, zz = abcd_to_xyz(abcd)
    for x, y, z in ((xx, yy, zz), (zz, xx, yy), (zz, yy, xx),):
        for sx in (x, -x):
            for sy in (y, -y):
                for sz in (z, -z):
                    u_set.add(xyz_to_u(sx, sy, sz))
    return u_set

"""
abcd_to_u_set is faster than abcd_to_u_set_uvw
r = repeat('abcd_to_u_set_uvw(abcd)', number=100, repeat=100, globals=globals()); np.median(r); np.std(r)
np.float64(0.018145770533010364)
np.float64(0.001256843665000622)
r = repeat('abcd_to_u_set(abcd)', number=100, repeat=100, globals=globals()); np.median(r); np.std(r)
np.float64(0.009961916599422693)
np.float64(0.001336732926780821)
"""

def known_to_bounds(known_limit: int=1000, 
                    m_limit: int=int(1e6), n_limit: int=int(1e7), ) -> tuple:
    """Search known solutions to estimate bound on u, v.
    
    Each u consists of a postive m=4k and a signed odd with with abs n.
    For each u in the set of 12 u for a solution,
    select the two u with smallest m and the 2 with smallest n.
    They will probably overlap.
    Keep the two n that are paired with the smallest m, 
    and the two m that are paired with the smallest n.
    The (m, n) bounds for u are the maxima of the first values,
    and the (m, n) bounds for v are the maxima of the second values.
    """
    um = [] # m paired with smallest n
    un = [] # n paired with smallest m
    vm = [] # m paired with 2nd smallest n
    vn = [] # n paired with 2nd smallest m
    for inx_known, val in enumerate(known.values(), start=1):
        abcd = val['abcd']
        u_set = abcd_to_u_set(abcd)
        mn = [0] * 12 # (m, n) for each u to be sorted
        nm = [0] * 12 # (n, m) for each u to be sorted
        for inx_u, u in enumerate(u_set):
            u = abs(u)
            m, n = numerator(u), denominator(u)
            if m % 2 != 0: m, n = n, m
            mn[inx_u] = m, n
            nm[inx_u] = n, m
        mn.sort()
        nm.sort()
        if nm[0][1] < m_limit:
            um.append(nm[0][1])
            if nm[1][1] < m_limit:
                vm.append(nm[1][1])
            else:
                print(f'Skipping vm {nm[1][1]} in #{inx_known}')
        else:
            print(f'Skipping um {nm[0][1]} and vm in #{inx_known}')

        if mn[0][1] < n_limit:
            un.append(mn[0][1])
            if mn[1][1] < n_limit:
                vn.append(mn[1][1])
            else:  
                print(f'Skipping vn {mn[1][1]} in #{inx_known}')
        else:  
            print(f'Skipping un {mn[0][1]} and vn in #{inx_known}')
        if len(um) >= known_limit:
            break
    return (f'{max(um):_}, {max(un):_}'), (f'{max(vm):_}, {max(vn):_}')
"""
known_to_bounds()
Skipping vm 19835764 in #10
Skipping un 11846053 and vn in #10
Skipping vm 12642040 in #54
Skipping vn 57878913 in #54
Skipping vm 12642040 in #56
Skipping vn 57878913 in #56
Skipping vm 10001951064 in #57
Skipping vn 8204718073 in #57
Skipping vm 29393447736 in #58
Skipping vn 14486729065 in #58
Skipping vm 4669000304 in #65
Skipping vn 944254963 in #65
Skipping vm 6985268 in #67
Skipping vn 92654145 in #67
Skipping vm 22869016 in #69
Skipping vm 116419537680 in #72
Skipping vn 3286280601165 in #72
Skipping vm 4669000304 in #73
Skipping vn 944254963 in #73
Skipping vm 6985268 in #74
Skipping vn 92654145 in #74
Skipping vm 323752247040 in #75
Skipping vn 634023165233 in #75
Skipping vm 1551044 in #76
Skipping vm 1856708 in #79
Skipping vm 29393447736 in #80
Skipping vn 14486729065 in #80
Skipping vm 61008600 in #84
Skipping vn 68433257 in #84
Skipping vm 4372152 in #85
Skipping vm 6037149201728 in #86
Skipping vn 3384184433553 in #86
Skipping vm 107014216 in #87
Skipping vn 210232185 in #87
Skipping vn 10490417 in #88
Skipping vm 116419537680 in #90
Skipping vn 22953456067 in #90
Skipping vm 541388136 in #92
Skipping vm 3885556 in #93
Skipping vn 23685689 in #93
Skipping vm 323752247040 in #94
Skipping vn 634023165233 in #94
('175_812, 93_017', '712_772, 4_037_701')
"""

def uv_to_Pk(u: Rational, v: Rational) -> tuple:
    """Calculate the parameters Pk from u, v
    """
    u2 = u * u
    u3 = u2 * u
    v2 = v * v
    v3 = v2 * v
    P0 =(2 + u2) * (2 + v2) * (12 - 8*u + 2*u2 - 8*v + 2*v2 + u2*v2)
    P1 = (-4 + 4*u + 2*v - u2*v) * (8*u -4*u2 + 4*v -8*u*v + 2*u2*v -4*v2 + 2*v3 + u2*v3)
    P2 = 2 * (-4 + 4*u + 2*v -u2*v) * (4 -2*u -4*v +u*v2) * (-u + v)
    P3 = (4 - 2*u - 4*v + u*v2) * (4*u - 4*u2 + 2*u3 + 8*v - 8*u*v - 4*v2 + 2*u*v2 + u3*v2)
    return (P0, P1, P2, P3)

def uvD_to_xyz(u: Rational, v: Rational, D: Rational) -> tuple:
    """Calculate (x, y, z) from u, v
    """
    P0, P1, P2, P3 = uv_to_Pk(u, v)
    denom = P0 + D*D
    u2 = u * u
    v2 = v * v
    x1 = (u2 - 2*u + 2) * (v2 - 2*v)
    y1 = 2*(u + v - 2) * (u + v - u*v)
    z1 = (v2 - 2*v + 2) * (u2 - 2*u)
    pair = [0]*2
    for inx, sD in enumerate((D, -D)):
        x = (P1 + x1 * sD) / denom
        y = (-(P1 + P2 + P3) + y1 * sD) / denom
        z = (P3 + z1 * sD) / denom
        pair[inx] = (x, y, z)
    return pair

def u_to_D_coeffs_int(um: int, un: int) -> tuple:
    """
    Returns int coefficients C0...C4 such that Ci / un^4 
    are the rational coefficients.
    """
    um2, un2 = um*um, un*un
    um3, un3 = um2*um, un2*un
    um4, un4 = um2*um2, un2*un2
 
    um_un3 = um * un3
    shared_um4_um3un = um4 - 4 * um3 * un
    
    c0 = 4 * (shared_um4_um3un + 8 * um_un3 - 12 * un4)
    c1 = 8 * (6 * um3 * un - um4 - 8 * um2 * un2 + 4 * um_un3 + 4 * un4)
    c2 = 16 * (3 * um2 * un2 - um3 * un - 4 * um_un3)
    c3 = 4 * (2 * um3 * un - um4 - 4 * um2 * un2 + 12 * um_un3 - 4 * un4)
    c4 = shared_um4_um3un - 8 * um_un3 + 4 * un4
    return [c0, c1, c2, c3, c4]

def v_to_D_int(vm: int, vn: int, Ci: tuple, un: int) -> None|Rational:
    """
    Evaluates the homogeneous polynomial V and checks if it's a square.
    V = C0*vn^4 + vm*(C1*vn^3 + vm*(C2*vn^2 + vm*(C3*vn + vm*C4)))
    """
    c0, c1, c2, c3, c4 = Ci
    vn2 = vn * vn
    vn3 = vn2 * vn
    vn4 = vn2 * vn2
    
    # Homogeneous Horner's Method
    V = c0 * vn4 + vm * (c1 * vn3 + vm * (c2 * vn2 + vm * (c3 * vn + vm * c4)))    
    if V < 0: return None
    # Bitwise filter: 75% of non-squares exit here in 1 cycle
    #if (V & 15) not in (0, 1, 4, 9): return None
    
    root, rem = Integer(V).sqrtrem()
    if rem == 0:
        # Reconstruct Rational D = sqrt(V) / (un^2 * vn^2)
        return QQ(root) / (un**2 * vn**2)
    return None

def search_generator_int(
        first_um: int, last_um: int, first_un: int, last_un: int,
        last_vm: int, last_vn: int):
    
    for um in range(first_um, last_um + 1, 4):
        print(f'um {um}')
        for un in range(first_un, last_un + 1, 2):
            if gcd(um, un) != 1: continue
            
            # (u_num, u_den) pairs for: u, -u, 1/u, -1/u
            u_variants = [(um, un), (-um, un), (un, um), (-un, um)]
            
            for u_num, u_den in u_variants:
                coeffs = u_to_D_coeffs_int(u_num, u_den)

                # Logic block 1: vn = um, vd > un
                vm = um
                for vn in range(un + 2, last_vn + 1, 2):
                    if gcd(vm, vn) != 1: continue
                    for v_num, v_den in [(vm, vn), (-vm, vn), (vn, vm), (-vn, vm)]:
                        D = v_to_D_int(v_num, v_den, coeffs, u_den)
                        if D is not None:
                            yield (QQ(u_num)/u_den, QQ(v_num)/v_den, D)
                
                # Logic block 2: vm > um
                for vm in range(um + 4, last_vm + 1, 4):
                    for vn in range(1, last_vn + 1, 2):
                        if gcd(vm, vn) != 1: continue
                        for v_num, v_den in [(vm, vn), (-vm, vn), (vn, vm), (-vn, vm)]:
                            D = v_to_D_int(v_num, v_den, coeffs, u_den)
                            if D is not None:
                                yield (QQ(u_num)/u_den, QQ(v_num)/v_den, D)

def find_solutions_int(first_um: int, last_um: int, first_un: int, last_un: int,
                   last_vm: int, last_vn: int,
                   max_d: int=int(1e27)) -> None | tuple:
    """Search for solutions with u over smaller range
    and v ranging from current u to larger range.
    The rationals u, v are constructed from 
        a postive multiple of 4 (m) and a signed odd value (n)
        The rational can be m/n, m/(-n), n/m, -n/m
    Breaks out of search if finds a new solution.
    """
    # Build search generator
    assert first_um%4 == last_um%4 == last_vm%4 == 0
    assert first_un%2 == last_un%2 == last_vn%2 == 1
    assert 4 <= first_um <= last_um < last_vm
    assert 1 <= first_un <= last_un < last_vn
    uv_est = int((last_um - first_um + 1) * (last_un - first_un + 1) 
              * last_vm * last_vn / (2.0 * 2.0)) # tends to be high

    # Look up index of known solution with denominator.
    d_to_known_inx = {val['abcd'][-1]: inx 
            for inx, val in enumerate(known.values(), start=1)}

    D_count = 0
    big_count = 0
    found_knowns = set() # indexes of known solutions found in search
    found_new = False
    start = time.time()

    for u, v, D in search_generator_int(first_um, last_um, first_un, last_un, last_vm, last_vn):
        if found_new: break
        print(f'{u}, {v} -> D {D}')
        D_count += 1
        pair = uvD_to_xyz(u, v, D)
        for xyz in pair:
            d = lcm([denominator(x) for x in xyz])
            if d >= max_d:
                big_count += 1
                print(f'\tbig {float(d):.4e}')
                continue
            if d in d_to_known_inx:
                inx = d_to_known_inx[d]
                found_knowns.add(inx)
                print(f'\tknown #{inx}: {xyz}')
                continue
            print(f'\n\nNEW ({u}, {v}) -> {xyz}')
            found_new = True
            break
    elapsed = time.time() - start
    print(f'uv_est {uv_est}, D_count {D_count}, big_count {big_count}'
          f'\nknowns {len(found_knowns)}: {found_knowns}')
    print(f'elapsed: {elapsed:.4f}s')
    if found_new:
        return u, v, xyz

"""
find_solutions_int(4, 20, 1, 9, 1000, 1041)
um 4
um 8
-5/8, -477/692 -> D 64335705/30647296
	known #7: (18796760/20615673, 2682440/20615673, -15365639/20615673)
	known #20: (2448718655/3393603777, -664793200/3393603777, -3134081336/3393603777)
um 12
um 16
um 20
-9/20, -1041/320 -> D 2126704839/40960000
	known #17: (1670617271/1679142729, 632671960/1679142729, -50237800/1679142729)
	known #1: (-95800/422481, -414560/422481, -217519/422481)
-9/20, 1000/47 -> D 495260031/441800
	known #1: (414560/422481, 95800/422481, 217519/422481)
	known #17: (-632671960/1679142729, -1670617271/1679142729, 50237800/1679142729)
uv_est 19909125, D_count 3, big_count 0
knowns 4: {17, 20, 1, 7}
elapsed: 81.8371s

find_solutions_int(4, 20, 1, 29, 1000, 1865)
um 4
um 8
-5/8, -1617/200 -> D 708019737/2560000
	known #14: (630662624/638523249, 275156240/638523249, 219076465/638523249)
	known #13: (-260052385/589845921, -582665296/589845921, -186668000/589845921)
-5/8, -477/692 -> D 64335705/30647296
	known #7: (18796760/20615673, 2682440/20615673, -15365639/20615673)
	known #20: (2448718655/3393603777, -664793200/3393603777, -3134081336/3393603777)
um 12
-29/12, 1865/132 -> D 4566509815/2509056
	known #34: (894416022327/2051764828361, -125777308440/2051764828361, 2032977944240/2051764828361)
	known #3: (-8332208/8707481, -5507880/8707481, -1705575/8707481)
um 16
um 20
-9/20, -1041/320 -> D 2126704839/40960000
	known #17: (1670617271/1679142729, 632671960/1679142729, -50237800/1679142729)
	known #1: (-95800/422481, -414560/422481, -217519/422481)
-9/20, -1425/412 -> D 3894577617/67897600
	known #23: (15355831360/15434547801, 5821981400/15434547801, -140976551/15434547801)
	known #2: (-673865/2813001, -2767624/2813001, -1390400/2813001)
-9/20, 1000/47 -> D 495260031/441800
	known #1: (414560/422481, 95800/422481, 217519/422481)
	known #17: (-632671960/1679142729, -1670617271/1679142729, 50237800/1679142729)
uv_est 114930625, D_count 6, big_count 0
knowns 10: {1, 34, 3, 2, 7, 13, 14, 17, 20, 23}
elapsed: 451.4831s

"""

def search_inner_loop_gen(um, first_un, last_un, last_vm, last_vn):
    print(f'um {um}')
    for un in range(first_un, last_un + 1, 2):
        print(f'un {un}')
        if gcd(um, un) != 1: continue
        
        # (u_num, u_den) pairs for: u, -u, 1/u, -1/u
        u_variants = [(um, un), (-um, un), (un, um), (-un, um)]
        
        for u_num, u_den in u_variants:
            coeffs = u_to_D_coeffs_int(u_num, u_den)

            # Logic block 1: vn = um, vd > un
            vm = um
            for vn in range(un + 2, last_vn + 1, 2):
                if gcd(vm, vn) != 1: continue
                for v_num, v_den in [(vm, vn), (-vm, vn), (vn, vm), (-vn, vm)]:
                    D = v_to_D_int(v_num, v_den, coeffs, u_den)
                    if D is not None:
                        yield (QQ(u_num)/u_den, QQ(v_num)/v_den, D)
            
            # Logic block 2: vm > um
            for vm in range(um + 4, last_vm + 1, 4):
                for vn in range(1, last_vn + 1, 2):
                    if gcd(vm, vn) != 1: continue
                    for v_num, v_den in [(vm, vn), (-vm, vn), (vn, vm), (-vn, vm)]:
                        D = v_to_D_int(v_num, v_den, coeffs, u_den)
                        if D is not None:
                            yield (QQ(u_num)/u_den, QQ(v_num)/v_den, D)

def persistent_worker(worker_id: int, task_queue, result_queue,
                      first_un, last_un, last_vm, last_vn):
    # This is the 'Local State' for this specific process
    while True:
        task_range = task_queue.get()
        if task_range is None: # The signal to stop
            result_queue.put(None)
            break
            
        start_u, end_u = task_range
        for um in range(start_u, end_u, -4):
            for r in search_inner_loop_gen(um, 
                    first_un, last_un, last_vm, last_vn):
                result_queue.put(r)

def chunk_generator(first_u: int, last_u: int, workers: int):
    """
    Ensures all workers get work by tapering down at the end.
    Allocates from end to beginning because larger u take longer.
    """
    current = last_u

    # Cap chunk size so one chunk doesn't eat the remainder.
    while current >= first_u:
        # Don't let a single chunk take more than remaining work per worker
        # and never exceed a reasonable max (like 50)
        remaining = (current - first_u) // 4 + 1
        max_chunk = max(1, min(remaining // workers, 50))
        yield (current, current - 4 * max_chunk)
        current -= 4 * max_chunk

def search_generator_p_int(first_um: int, last_um: int,
                first_un: int, last_un: int, last_vm: int, last_vn: int,
                workers: int=8):
    """
    Performs the search in parallel using a dynamic work pool.
    """
    task_queue = Queue()
    result_queue = Queue()

    processes = []
    for i in range(workers):
        p = Process(target=persistent_worker,
                    args=(i, task_queue, result_queue, # Pass both queues
                          first_un, last_un, last_vm, last_vn))
        p.start()
        processes.append(p)
        
    for chunk in chunk_generator(first_um, last_um, workers):
        task_queue.put(chunk)
        
    for _ in range(workers):
        task_queue.put(None)

    # Listen for results while workers are running
    finished_workers = 0
    while finished_workers < workers:
        result = result_queue.get() 
        if result is None:
            finished_workers += 1
        else:
            yield result # Send (u, v, D) back to the main function

    # Cleanup
    for p in processes:
        p.join()

def find_solutions_p_int(first_um: int, last_um: int, first_un: int, last_un: int,
                   last_vm: int, last_vn: int, workers: int=8,
                   max_d: int=int(1e27)) -> None | tuple:
    """Parallel search for solutions with u over smaller range
    and v ranging from current u to larger range.
    The rationals u, v are constructed from 
        a postive multiple of 4 (m) and a signed odd value (n)
        The rational can be m/n, m/(-n), n/m, -n/m
    Breaks out of search if finds a new solution.
    """
    # Build search generator
    assert first_um%4 == last_um%4 == last_vm%4 == 0
    assert first_un%2 == last_un%2 == last_vn%2 == 1
    assert 4 <= first_um <= last_um < last_vm
    assert 1 <= first_un <= last_un < last_vn
    uv_est = int((last_um - first_um + 1) * (last_un - first_un + 1) 
              * last_vm * last_vn / (2.0 * 2.0)) # tends to be high

    # Look up index of known solution with denominator.
    d_to_known_inx = {val['abcd'][-1]: inx 
            for inx, val in enumerate(known.values(), start=1)}

    D_count = 0
    big_count = 0
    found_knowns = set() # indexes of known solutions found in search
    found_new = False
    start = time.time()

    for u, v, D in search_generator_p_int(
            first_um, last_um, first_un, last_un, last_vm, last_vn, workers):
        print(f'{u}, {v} -> D {D}')
        D_count += 1
        pair = uvD_to_xyz(u, v, D)
        for xyz in pair:
            d = lcm([denominator(x) for x in xyz])
            if d >= max_d:
                big_count += 1
                print(f'\tbig {float(d):.4e}')
                continue
            if d in d_to_known_inx:
                inx = d_to_known_inx[d]
                found_knowns.add(inx)
                print(f'\tknown #{inx}: {xyz}')
                continue
            print(f'\n\nNEW ({u}, {v}) -> {xyz}')
            found_new = True
            break
        if found_new: break
    elapsed = time.time() - start
    print(f'uv_est {uv_est}, D_count {D_count}, big_count {big_count}'
          f'\nknowns {len(found_knowns)}: {found_knowns}')
    print(f'elapsed: {elapsed:.4f}s')
    if found_new:
        return u, v, xyz

def known_to_big_v(lim_d_known: int=int(1e27), lim_d: int=int(1e40)) -> tuple:
    """Use known solutions to generate sets of v,
    and use them to generate new solutions below lim_d,
    and use them to generate a big set of v.
    If any of the new solutions are below lim_d_known
    and not in set of known, report new solution.
    
    Return big_v_set, known_denoms  
    """    
    known_denoms = set(val['abcd'][-1] for val in known.values())
    denoms = set()
    denoms.update(known_denoms)
    big_v_set = set()
    for val in known.values():
        # use known solutions to generate set of v,
        abcd = val['abcd']
        u_set = abcd_to_u_set(abcd)
        assert 12 == len(u_set)
        big_v_set.update(u_set)

        # use v_set to generate solutions
        v_n_d_pairs = [(numerator(v), denominator(v)) for v in u_set]
        for i, (u_num, u_den) in enumerate(v_n_d_pairs[:-1]):
            coeffs = u_to_D_coeffs_int(u_num, u_den)
            for v_num, v_den in v_n_d_pairs[i+1:]:
                D = v_to_D_int(v_num, v_den, coeffs, u_den)
                if D is None: continue
                pair = uvD_to_xyz(QQ(u_num)/u_den, QQ(v_num)/v_den, D)
                assert pair is not None
                for xyz in pair:
                    d = lcm([denominator(x) for x in xyz])
                    if d >= lim_d: continue
                    if d in denoms: continue
                    denoms.add(d)
                    abc = sorted(abs(x) * d for x in xyz)
                    abc_d = abc + [d]
                    if d < lim_d_known:
                        assert d in known_denoms, f'New solution: {QQ(u_num)/u_den}, {QQ(v_num)/v_den} -> {abc_d}'
                    v_set = abcd_to_u_set(abc_d)
                    assert 12 == len(v_set)
                    big_v_set.update(v_set)
    # Restrict v_set to numerator and denominator < lim_d
    v_set = set()
    for v in big_v_set:
        if max(abs(numerator(v)), abs(denominator(v))) < lim_d:
            v_set.add(v)
    return big_v_set, known_denoms

"""
% python -m elliptical.brute known_v_to_brute_u 
1003 without big
1723 with big
"""
def known_v_to_brute_u(min_umn: int=1, max_umn=10, 
                    lim_d_known: int=int(1e27), lim_d: int=int(1e40)) -> None|tuple:
    """Use known solutions to generate a big set of v,
    Select the v with either numerator or denominator above max_unm,
    and brute force search on these v with u in given umn range.
    Break out of search if find a new solution below max_d_known.
    """
    big_v_set, known_denoms = known_to_big_v(lim_d_known, lim_d)
    use_v_set = set()
    for v in big_v_set:
        v_num = numerator(v)
        v_den = denominator(v)
        assert 0 < v_den
        v_max = max(abs(v_num), v_den)
        if  max_umn < v_max < lim_d:
            use_v_set.add((v_num, v_den))
    print(f'Using {len(use_v_set)} of {len(big_v_set)} v')

    assert 1 <= min_umn <= max_umn 
    
    #first_um = 4 + min_umn - (min_umn & 0b11) if min_umn & 0b11 else min_umn
    min_un = min_umn | 0b1
    last_um = last_un = max_umn

    denoms = set()
    denoms.update(known_denoms)
    D_count = 0
    big_count = 0
    known_count = 0
    for um in range(4, last_um + 1, 4):
        print(f'um {um}')
        first_un = min_un if um < min_un else 1
        for un in range(first_un, last_un + 1, 2):
            if gcd(um, un) != 1: continue
            #print(f'un {un}')
            u_variants = [(um, un), (-um, un), (un, um), (-un, um)]            
            for u_num, u_den in u_variants:
                coeffs = u_to_D_coeffs_int(u_num, u_den)
                for v_num, v_den in use_v_set:
                    D = v_to_D_int(v_num, v_den, coeffs, u_den)
                    if D is None: continue
                    D_count += 1
                    u = QQ(u_num)/u_den
                    v = QQ(v_num)/v_den
                    pair = uvD_to_xyz(u, v, D)
                    for xyz in pair:
                        d = lcm([denominator(x) for x in xyz])
                        if d >= lim_d_known:
                            big_count += 1
                            continue
                        if d in known_denoms:
                            known_count += 1
                            continue
                        if d in denoms: continue
                        denoms.add(d)
                        # FOUND NEW !!!!!!!
                        abc = sorted(abs(x) * d for x in xyz)
                        abc_d = abc + [d]
                        assert d in known_denoms, f'New solution: {QQ(u_num)/u_den}, {QQ(v_num)/v_den} -> {abc_d}'
    print(f'V count {len(use_v_set)}, D_count {D_count}, big_count {big_count}, known_count {known_count}')

"""
V count 1317, D_count 36, big_count 32, known_count 40
python -m elliptical.brute known_v_to_brute_u 1 10  1.78s user 0.27s system 99% cpu 2.060 total

V count 1312, D_count 138, big_count 128, known_count 148
python -m elliptical.brute known_v_to_brute_u 1 100  8.73s user 0.29s system 99% cpu 9.040 total

V count 1312, D_count 102, big_count 96, known_count 108
python -m elliptical.brute known_v_to_brute_u 11 100  8.62s user 0.29s system 99% cpu 8.919 total

V count 1299, D_count 256, big_count 228, known_count 284
python -m elliptical.brute known_v_to_brute_u 1 1000  691.66s user 1.96s system 99% cpu 11:34.64 total

V count 1288, D_count 101, big_count 112, known_count 90
python -m elliptical.brute known_v_to_brute_u 1001 2000  2057.17s user 4.57s system 99% cpu 34:23.39 total
"""

def search_v_to_brute_u_gen(um: int, min_un: int, 
                            last_um: int, last_un: int,
                            use_v_set: set):
    print(f'um {um}')
    first_un = min_un if um < min_un else 1
    for un in range(first_un, last_un + 1, 2):
        if gcd(um, un) != 1: continue
        u_variants = [(um, un), (-um, un), (un, um), (-un, um)]            
        for u_num, u_den in u_variants:
            coeffs = u_to_D_coeffs_int(u_num, u_den)
            for v_num, v_den in use_v_set:
                D = v_to_D_int(v_num, v_den, coeffs, u_den)
                if D is not None:
                        yield (QQ(u_num)/u_den, QQ(v_num)/v_den, D)

def persistent_worker_v_to_brute_u(worker_id: int, task_queue, result_queue,
                min_un: int, last_um: int, last_un: int, use_v_set: set):
    # This is the 'Local State' for this specific process
    while True:
        task_range = task_queue.get()
        if task_range is None: # The signal to stop
            result_queue.put(None)
            break
            
        start_u, end_u = task_range
        for um in range(start_u, end_u, -4):
            for r in search_v_to_brute_u_gen(um, min_un, 
                            last_um, last_un, use_v_set):
                result_queue.put(r)

def search_generator_v_to_brute_u(min_un: int, last_um: int, last_un: int, 
                                  use_v_set: set, workers: int=8):
    """ Parallel search generator for known_v_to_brute_u_p .
    """
    task_queue = Queue()
    result_queue = Queue()

    processes = []
    for i in range(workers):
        p = Process(target=persistent_worker_v_to_brute_u,
                    args=(i, task_queue, result_queue, # Pass both queues
                          min_un, last_um, last_un, use_v_set))
        p.start()
        processes.append(p)
        
    #for chunk in chunk_generator(4, min_un, workers): # patch to run only remainder quadrant
    for chunk in chunk_generator(4, last_um, workers):
        task_queue.put(chunk)
        
    for _ in range(workers):
        task_queue.put(None)

    # Listen for results while workers are running
    finished_workers = 0
    while finished_workers < workers:
        result = result_queue.get() 
        if result is None:
            finished_workers += 1
        else:
            yield result # Send (u, v, D) back to the main function

    # Cleanup
    for p in processes:
        p.join()

def known_v_to_brute_u_p(min_umn: int=1, max_umn=10, 
                    lim_d_known: int=int(1e27), lim_d: int=int(1e40)) -> None|tuple:
    """ Parallel version of known_v_to_brute_u .
    Use known solutions to generate a big set of v,
    Select the v with either numerator or denominator above max_unm,
    and brute force search on these v with u in given umn range.
    Break out of search if find a new solution below max_d_known.
    """
    big_v_set, known_denoms = known_to_big_v(lim_d_known, lim_d)
    use_v_set = set()
    for v in big_v_set:
        v_num = numerator(v)
        v_den = denominator(v)
        assert 0 < v_den
        v_max = max(abs(v_num), v_den)
        if  max_umn < v_max < lim_d:
            use_v_set.add((v_num, v_den))
    print(f'Using {len(use_v_set)} of {len(big_v_set)} v')

    assert 1 <= min_umn <= max_umn 
    min_un = min_umn | 0b1
    last_um = last_un = max_umn

    D_count = 0
    big_count = 0
    known_count = 0
    for u, v, D in search_generator_v_to_brute_u(
            min_un, last_um, last_un, use_v_set):
        D_count += 1
        pair = uvD_to_xyz(u, v, D)
        for xyz in pair:
            d = lcm([denominator(x) for x in xyz])
            if d >= lim_d_known:
                big_count += 1
                continue
            if d in known_denoms:
                known_count += 1
                continue
            if d in denoms: continue
            denoms.add(d)
            # FOUND NEW !!!!!!!
            abc = sorted(abs(x) * d for x in xyz)
            abc_d = abc + [d]
            assert d in known_denoms, f'New solution: {QQ(u_num)/u_den}, {QQ(v_num)/v_den} -> {abc_d}'
    print(f'V count {len(use_v_set)}, D_count {D_count}, big_count {big_count}, known_count {known_count}')

"""
V count 1299, D_count 256, big_count 228, known_count 284
python -m elliptical.brute known_v_to_brute_u_p 1 1000  749.37s user 4.07s system 764% cpu 1:38.54 total
"""

def uv_to_solutions_explore(u_e, s_u_o, v_e, s_v_o):
    """Explore which orderings of parts of u and v generate unique solutions.
    """
    for u_num, u_den in ((u_e, s_u_o), (s_u_o, u_e)):
        coeffs = u_to_D_coeffs_int(u_num, u_den)
        for v_num, v_den in ((v_e, s_v_o), (s_v_o, v_e)):
            D = v_to_D_int(v_num, v_den, coeffs, u_den)
            if D is None:
                print(u_num, u_den, v_num, v_den, 'FAIL')
                continue
            u = QQ(u_num)/u_den
            v = QQ(v_num)/v_den
            pair = uvD_to_xyz(u, v, D)
            for xyz in pair:
                d = lcm([denominator(x) for x in xyz])
                print(u_num, u_den, v_num, v_den, d)

"""
Only one ordering produces two solutions. And both can be negative.
>>> uv_to_solutions_explore(20, -9, 320, -1041)
20 -9 320 -1041 FAIL
20 -9 -1041 320 FAIL
-9 20 320 -1041 FAIL
-9 20 -1041 320 1679142729
-9 20 -1041 320 422481

And other signs fail including both positive
>>> uv_to_solutions_explore(20, -9, 320, 1041)
>>> uv_to_solutions_explore(20, 9, 320, 1041)
>>> uv_to_solutions_explore(20, 9, 320, -1041)
"""

def make_uv_table(max_umn: int, max_vmn: int,
                  lim_d_known: int=int(1e27), lim_d: int=int(1e40)) -> None:
    """Make a table of all the extended known u, v and the solutions they generate
    """
    big_v_set, known_denoms = known_to_big_v(lim_d_known, lim_d)
    print(f'big_v_set has {len(big_v_set)}')

    # Separate the even and odd parts of the v
    evens = set() # set of abs values of even numerators, denominators
    odds = set() # set of abs values of odd numerators, denominators
    for v in big_v_set:
        for n in (numerator(v), denominator(v)):
            n = abs(n)
            if n % 2 == 0: evens.add(n)
            else: odds.add(n)
    print(f'even parts {len(evens)} up to {max(evens):.2e}')
    print(f'odd parts {len(odds)} up to {max(odds):.2e}')

    # Generate each possible rational from the evens and odds and their solutions
    evens = sorted(evens)
    u_evens_max_inx = searchsorted(evens, max_umn, 'left')
    v_evens_max_inx = searchsorted(evens, max_vmn, 'left')
    odds = sorted(odds)
    u_odds_max_inx = searchsorted(odds, max_umn, 'left')
    v_odds_max_inx = searchsorted(odds, max_vmn, 'left')
    count = 0
    D_count = 0
    with open('solutions_uv.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['u_num', 'u_den', 'v_num', 'v_den', 'd', 'c', 'b', 'a'])
        for u_inx_e, u_e in enumerate(evens[: u_evens_max_inx]):
            for u_o in odds[:u_odds_max_inx]:
                if gcd(u_e, u_o) != 1: continue
                for u_num, u_den in ((-u_o, u_e), (u_e, -u_o), (u_o, u_e), (u_e, u_o)):
                    coeffs = u_to_D_coeffs_int(u_num, u_den)
                    for v_e in evens[u_inx_e + 1: v_evens_max_inx]:
                        for v_o in odds[: v_evens_max_inx]:
                            if gcd(v_e, v_o) != 1: continue
                            for v_num, v_den in ((-v_o, v_e), (v_e, -v_o), (v_o, v_e), (v_e, v_o)):
                                count +=1
                                D = v_to_D_int(v_num, v_den, coeffs, u_den)
                                if D is None: continue
                                D_count += 1
                                u = QQ(u_num)/u_den
                                v = QQ(v_num)/v_den
                                pair = uvD_to_xyz(u, v, D)
                                for xyz in pair:
                                    d = int(lcm([denominator(x) for x in xyz]))
                                    c, b, a = sorted([abs(int(x * d)) for x in xyz], reverse=True)
                                    row = [u_num, u_den, v_num, v_den, d, c, b, a]
                                    writer.writerow(row)
                                    if d < lim_d_known:
                                        assert d in known_denoms, f'New solution: {row}'
    print(count, D_count)
    
""" 
% python -m elliptical.brute make_uv_table 1_000 1_000
<function make_uv_table at 0x174075800>
big_v_set has 1723
even parts 1720 up to 4.33e+77
odd parts 1722 up to 4.29e+76
2231008 6
cpu 3.714 total

u_num,u_den,v_num,v_den,d,c,b,a
201,4,136,-133,156646737,146627384,108644015,27450160
201,4,136,-133,16003017,14173720,12552200,4479031
201,4,-1005,568,16003017,14173720,12552200,4479031
201,4,-1005,568,156646737,146627384,108644015,27450160
-5,8,-477,692,20615673,18796760,15365639,2682440
-5,8,-477,692,3393603777,3134081336,2448718655,664793200
-9,20,-1041,320,1679142729,1670617271,632671960,50237800
-9,20,-1041,320,422481,414560,217519,95800
-93,80,400,-37,12197457,11289040,8282543,5870000
-93,80,400,-37,1787882337,1662997663,1237796960,686398000
136,-133,-1005,568,156646737,146627384,108644015,27450160
136,-133,-1005,568,16003017,14173720,12552200,4479031

% time python -m elliptical.brute make_uv_table 1_000_000 12_000_000
900173248 132
cpu 9:59.07 total

% time python -m elliptical.brute make_uv_table 12_000_000 100_000_000
2840930576 191
cpu 32:04.23 total
"""


if __name__ == "__main__":
    import sys
    command = eval(sys.argv[1])
    print(command)
    args = list(map(int, sys.argv[2:]))
    result = command(*args)
    print(result)
