"""
Sequence utilities for LENSTools.
"""

def iupac_conversion(base):
    """
    Convert a base to A, C, T, G space.

    Args:
        base (str): Input base

    Returns:
        emission (str): Output base
    """
    emission = []
    conversion = {'R': ['A', 'G'],
                  'Y': ['C', 'T'],
                  'S': ['G', 'C'],
                  'W': ['A', 'T'],
                  'K': ['G', 'T'],
                  'M': ['A', 'C']}
    if base in conversion.keys():
        emission = conversion[base]
    return emission
