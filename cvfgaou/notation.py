"""Helper constants and functions to maintain consistent entity notation."""

import pandas as pd

# Mapping of AA names
AA_NAMES = {
    'A': 'Ala',
    'R': 'Arg',
    'N': 'Asn',
    'D': 'Asp',
    'C': 'Cys',
    'Q': 'Gln',
    'E': 'Glu',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'L': 'Leu',
    'K': 'Lys',
    'M': 'Met',
    'F': 'Phe',
    'P': 'Pro',
    'S': 'Ser',
    'T': 'Thr',
    'W': 'Trp',
    'Y': 'Tyr',
    'V': 'Val',
    'U': 'Sec',
    'O': 'Pyl'
}

# Unicode character codes
GEQ_CHAR = '\u2265'
LEQ_CHAR = '\u2264'

# Strength of evidence naming/reference table
# Stores correspondence of point magnitude to strength phrase
SOE_TABLE = (
    (8, 'very strong'),
    (4, 'strong'),
    (3, 'moderate+'),
    (2, 'moderate'),
    (1, 'supporting')
)

# Direction of evidence naming/reference table
# Stores directional information:
# inequality, sign, criterion, label, comparator, reverse_comparator
DOE_TABLE = (
    (GEQ_CHAR, '+', 'PP3', 'Pathogenic', pd.Series.ge, pd.Series.le),
    (LEQ_CHAR, '-', 'BP4', 'Benign', pd.Series.le, pd.Series.ge)
)

evidence_phrases_to_points = {
    f'{criterion} {strength.title()}'.lower(): points * sign
    for criterion, sign in (('PP3', 1), ('BP4', -1))
    for points, strength in SOE_TABLE
}

def evidence_meets_strength(evidence, reference, strict=False):
    """Compare an evidence label to a reference point

    `evidence` and `reference` are reresentations of directed evidence
    strength; these can be strings such as "PP3 very strong", "+8" or
    integers.
    Returns true if evidence meets the direction of reference and is at
    least as strong, i.e. (+8, +4) and (-8, -4) should return True, and
    (+8, -4) should return False.
    If `strict` is True, raises a ValueError if `evidence` cannot be
    parsed, but if `strict` is False, failure to parse `evidence` will
    simply return False for all valid values of `reference`.
    """

    if isinstance(reference, str):
        ref_lookup = evidence_phrases_to_points.get(reference.replace('_', ' ').lower())
        if ref_lookup is None:
            try:
                reference = int(reference)
            except:
                raise ValueError(f'Could not parse {reference} as an evidence strength.')
        else:
            reference = ref_lookup
    
    if reference == 0: raise ValueError('Reference evidence strength cannot be 0')

    if isinstance(evidence, str):
        ev_lookup = evidence_phrases_to_points.get(evidence.replace('_', ' ').lower())
        if ev_lookup is None:
            try:
                evidence = int(evidence)
            except:
                if strict:
                    raise ValueError(f'Could not parse {evidence} as an evidence strength')
                else:
                    return False
    
    return evidence >= reference if reference > 0 else evidence <= reference


def var_str(contig, pos, ref, alt):
    """Compose the variant ID string contig:pos:ref>alt"""
    return f'{contig}:{pos}:{ref}>{alt}'


def var_col(df, contig_col, pos_col, ref_col, alt_col):
    """Create a variant ID column given contig, pos, ref, alt columns in a dataframe"""
    return df.apply(
        lambda row: var_str(
            row[contig_col],
            row[pos_col],
            row[ref_col],
            row[alt_col]
        ),
        axis='columns'
    )