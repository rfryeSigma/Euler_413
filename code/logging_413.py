# Common logging class.
import argparse
from enum import IntFlag, auto

class V(IntFlag):
    # Flags for selecting logging messages
    # Defining flags using auto() assigns powers of 2 (1, 2, 4, 8, ...)
    NONE   = 0x00000
    DEBUG  = auto() # 0x00001
    SEARCH = 0x0010
    OUTER  = auto() # 0x00020
    INNER  = auto() # 0x00040
    PROCESS= auto() # 0x00080
    FACTOR = 0x00100
    F_COMMON=auto() # 0x00200
    F_P_1  = auto() # 0x00400
    F_GCD  = auto() # 0x00800
    F_TRIAL= auto() # 0x01000
    F_RHO  = auto() # 0x02000
    F_DIAG = auto() # 0x04000
    SOLVE  = 0x10000
    S_FACT = auto() # 0x20000
    S_PART = auto() # 0x40000
    S_CUBIC= auto() # 0x80000
    ALL = (2**(S_CUBIC.bit_length())) - 1
    # END Flags for selecting logging messages
    # Flags for selecting actions
    RHO1 = auto() # 0x100000 # call pollard_rho_one
    
    @staticmethod
    def log(flags, target, message, print=print):
        """Prints message only if target is set within flags."""
        if flags and flags & target: print(message)


def parse_flags(flag_string):
    """Converts command line string into the combined IntFlag value."""
    try:
        # Split by '|' and sum the attributes found in V
        result = V.NONE
        for part in flag_string.split('|'):
            result |= V[part.upper().strip()]
        return result
    except KeyError as e:
        raise argparse.ArgumentTypeError(f"Invalid flag: {e.args[0]}. Valid: {list(V.__members__.keys())}")
