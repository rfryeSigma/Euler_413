# Utility functions used during the 413 process

from math import gcd

def load_tree(level: int, folder: str) -> list:
    """
    Load all or a single level of a prime product tree.
    The 0th level is the leaves, the primes themselves.
    Intermediate levels are pairwise products of previous level
    The last level is a list containing the product of all primes.
    
    Args:
        level (int): level to load;
            if -1, load root;
            if -2 load all levels
        folder (str): name of folder containing level files
    
    Returns: 
        list of products at selected level or list of lists.
    """
    def load_level(level_idx: int) -> list:
        """Loads a specific level into memory."""
        num_nodes = 2** (num_levels - level_idx - 1)
        nodes = [[]] * num_nodes
        path = f"{folder}/level_{level_idx}.bin"
        with open(path, "rb") as f:
            for node_idx in range(num_nodes):
                len_bytes = f.read(4)
                if not len_bytes:
                    return nodes[:node_idx]
                    break
                length = int.from_bytes(len_bytes, 'little')
                nodes[node_idx] = (int.from_bytes(f.read(length), 'little'))
            else:
                assert node_idx + 1 == num_nodes
            return nodes

    def get_root():
        """Returns the single product at the top of the tree."""
        return load_level(num_levels - 1)[0]

    with open(f"{folder}/meta.txt", "r") as f:
        num_levels = int(f.read())
        #print(f'{folder} has {num_levels} levels')

        if 0 <= level:
            nodes = load_level(level)
            return nodes
        
        if -1 == level:
            return get_root()
        
        all_nodes = [[]] * num_levels
        for level in range(num_levels):
            all_nodes[level] = load_level(level)
        return all_nodes

def factor_tree_gen(g: int, tree: list, level: int | None = None, index: int = 0):
    """
    Yield prime factors of g by recursively descending prime product tree.
    """
    if level is None:
        level = len(tree) - 1

    # Base case: we reached the leaves
    if level == 0:
        p = tree[0][index]
        if g % p == 0: yield p
        return

    # Check the "Left" child
    left_index = index * 2
    left_prod = tree[level-1][left_index]
    g_left = gcd(g, left_prod)
    if g_left > 1:
        yield from factor_tree_gen(g_left, tree, level-1, left_index)

    # Check the "Right" child
    right_index = index * 2 + 1
    if right_index < len(tree[level-1]):
        right_prod = tree[level-1][right_index]
        g_right = gcd(g, right_prod)
        if g_right > 1:
            yield from factor_tree_gen(g_right, tree, level-1, right_index)
