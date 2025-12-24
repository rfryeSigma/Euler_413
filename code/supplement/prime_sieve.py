# Hand code single buffer segmented sieve of Eratosthenes.

import math
import os
import time
import argparse

# --- Global Constants (from prompt) ---
N_PRIMES = int(1e6) # Generate at least this many primes
BUFFER_OVERHEAD = int(2**16) # must be >= math.ceil(math.ln(max_prime)**2) 

TOTAL_BUFFER_SIZE = N_PRIMES + BUFFER_OVERHEAD
initial_primes = (2, 3,)  # kept minimal; other primes found by sieve
prime_count = len(initial_primes)

# --- Global Buffer Declaration (Constraint) ---
# prime_buffer stores discovered primes in its low indices; the segmented
# bytearray `segment_buf` is reused for marking odd candidates in each pass.
prime_buffer = [0] * TOTAL_BUFFER_SIZE
segment_buf = bytearray(BUFFER_OVERHEAD)

def simple_sieve():
    """
    Uses primes currently in prime_buffer to sieve the next segment.
    """
    # Declare mutated globals and localize for speed
    global prime_count
    pb = prime_buffer
    seg = segment_buf

    high_prime = pb[prime_count - 1]
    next_limit = high_prime * high_prime
    # number of odd candidates between high_prime and next_limit (exclusive)
    segment_len = (next_limit - high_prime) // 2 - 1
    # don't exceed reserved buffer capacity or the reusable segment buffer
    max_new = TOTAL_BUFFER_SIZE - prime_count
    seg_len_cap = len(seg)
    if segment_len > max_new:
        segment_len = max_new
    if segment_len > seg_len_cap:
        segment_len = seg_len_cap
    if segment_len <= 0:
        return

    # initialize segment bytes to 1 (candidate primes)
    seg_view = memoryview(seg)
    seg_view[:segment_len] = b'\x01' * segment_len

    # Pre-filter multiples of small primes (3 and 5) to reduce marking
    for small_p in (3, 5):
        first_multiple = ((high_prime + small_p) // small_p) * small_p
        # start marking at least from p*p to avoid clearing the prime itself
        if first_multiple < small_p * small_p:
            first_multiple = small_p * small_p
        if first_multiple % 2 == 0:
            first_multiple += small_p
        start = (first_multiple - high_prime - 2) // 2
        if 0 <= start < segment_len:
            num = (segment_len - 1 - start) // small_p + 1
            if num > 0:
                seg[start:segment_len:small_p] = b'\x00' * num

    # sieve using previously found primes (skip 2, and skip 3/5 marking)
    for inx_p in range(1, prime_count):
        p = pb[inx_p]
        if p == 3 or p == 5:
            continue
        # Find first odd multiple of p greater than high_prime and >= p*p
        first_multiple = ((high_prime + p) // p) * p
        if first_multiple < p * p:
            first_multiple = p * p
        if first_multiple % 2 == 0:
            first_multiple += p
        start = (first_multiple - high_prime - 2) // 2
        if start >= segment_len:
            continue
        step = p
        # mark composites in the reused bytearray using slice assignment
        if start < segment_len:
            num = (segment_len - 1 - start) // step + 1
            if num > 0:
                seg[start:segment_len:step] = b'\x00' * num

    # compact new primes into prime_buffer
    inx_b = prime_count
    base = high_prime
    for i in range(segment_len):
        if seg[i]:
            val = base + 2 * (i + 1)
            pb[inx_b] = val
            inx_b += 1
            # stop exactly at requested N_PRIMES to avoid overshoot
            if inx_b >= N_PRIMES:
                prime_count = N_PRIMES
                return
            if inx_b >= TOTAL_BUFFER_SIZE:  # safety
                break
    prime_count = inx_b

def run_sieve(verbose: bool = False):
    """
    Fill prime_buffer by sieving segments until at least N_PRIMES primes are generated.
    """
    if verbose:
        print(f"Generate at least {N_PRIMES} primes with buffer overhead {BUFFER_OVERHEAD}...")
    prime_buffer[:prime_count] = initial_primes
    start = time.time()
    seg_count = 0
    log_every = 10
    while prime_count < N_PRIMES:
        simple_sieve()
        seg_count += 1
        # log progress periodically
        if verbose and (seg_count % log_every == 0 or prime_count >= N_PRIMES):
            elapsed = time.time() - start
            pct = (prime_count / N_PRIMES) * 100 if N_PRIMES else 0
            eta = None
            if prime_count > 0:
                try:
                    eta = elapsed * (N_PRIMES / prime_count - 1)
                except ZeroDivisionError:
                    eta = None
            eta_s = f"{eta:.1f}s" if eta is not None else "unknown"
            print(f"[sieve] segments={seg_count} primes={prime_count}/{N_PRIMES} "
                  f"({pct:.2f}%) elapsed={elapsed:.1f}s ETA={eta_s}")
    return prime_count

def select_1mod8_primes(n_primes=2**17, start_inx=0) -> list:
    """
    Select primes from prime_buffer such that p%8 == 1.

    
    :param n_primes: number of primes to select; should be power of 2 if tree.
    :param start_inx: where to start in prime buffer
    """
    primes = [0] * n_primes
    buf_inx = start_inx
    for p_inx in range(n_primes):
        while prime_buffer[buf_inx] % 8 != 1:
            buf_inx += 1
        primes[p_inx] = prime_buffer[buf_inx]
        buf_inx += 1
    return primes

def build_tree(tree_primes: list) -> list:
    """
    Build a product tree for a list of primes
    
    :param tree_primes: list of primes to tree; best if length power of 2.
    """
    # Build Tree (Pass-through for odd levels)
    tree = [tree_primes]
    while len(tree[-1]) > 1:
        curr = tree[-1]
        next_lvl = [curr[i] * curr[i+1] for i in range(0, len(curr)-1, 2)]
        if len(curr) % 2 == 1: next_lvl.append(curr[-1])
        tree.append(next_lvl)
    return tree

def save_prime_tree(primes, folder="prime_tree"):
    """
    Construct a product tree for the selected primes.
    
    :param primes: selection of primes to save. Ideal length is a power of two.
    :param folder: name of folder to hold the prime tree.
    """
    if len(primes).bit_count() != 1:
        print(f'WARNING: number of primes {len(primes)} is not a power of 2')
    tree = build_tree(primes)

    if not os.path.exists(folder):
        os.makedirs(folder)

    # Save leaf primes (level 0) and internal levels
    for i in range(len(tree)):
        with open(f"{folder}/level_{i}.bin", "wb") as f:
            for node in tree[i]:
                node_bytes = node.to_bytes((node.bit_length() + 7) // 8, 'little')
                f.write(len(node_bytes).to_bytes(4, 'little'))
                f.write(node_bytes)
    
    # Save a metadata file so we know how many levels there are
    with open(f"{folder}/meta.txt", "w") as f:
        f.write(str(len(tree)))

def load_prime_tree(level=None, folder="prime_tree") -> list:
    """
    Load all or a single level of a prime product tree.
    The 0th level is the leaves, the primes themselves.
    Intermediate levels are pairwise products of previous level
    The last level is a list containing the product of all primes.
    
    :param level: level to load; if None, load all levels
    :param folder: name of folder containing level files
    :return: list of products at selected level or list of lists.
    :rtype: list
    """
    with open(f"{folder}/meta.txt", "r") as f:
        num_levels = int(f.read())
        print(f'{folder} has {num_levels} levels')

    def load_level(level_idx):
        """Loads a specific level into memory only when needed."""
        num_nodes = 2** (num_levels - level_idx - 1)
        nodes = [[]] * num_nodes
        path = f"{folder}/level_{level_idx}.bin"
        with open(path, "rb") as f:
            for node_idx in range(num_nodes):
                len_bytes = f.read(4)
                print(node_idx, len_bytes)
                if not len_bytes: 
                    print(f'break after build {node_idx} nodes out of {num_nodes}')
                    break
                length = int.from_bytes(len_bytes, 'little')
                nodes[node_idx] = (int.from_bytes(f.read(length), 'little'))
        return nodes

    def get_root():
        """Returns the single product at the top of the tree."""
        return load_level(num_levels - 1)[0]

    if level is not None:
        nodes = load_level(level)
        return nodes
    
    all_nodes = [[]] * num_levels
    for level in range(num_levels):
        all_nodes[level] = load_level[level]
    return all_nodes


# --- Main Section ---
def main(argv=None):
    parser = argparse.ArgumentParser(description='Run prime sieve')
    parser.add_argument('-v', '--verbose', action='store_true', help='enable segment logging')
    args = parser.parse_args(argv)
    count = run_sieve(verbose=args.verbose)
    assert count == prime_count

    # Display results and verify constraints
    print_len = 10
    print(f"Generated {prime_count} primes")
    print(f"\t{prime_buffer[:print_len]} ...")
    print(f"\t{prime_buffer[(N_PRIMES - print_len):N_PRIMES]} ...")
    print(f"\t{prime_buffer[prime_count -1]}")


if __name__ == "__main__":
    import sys
    raise SystemExit(main(sys.argv[1:]))
