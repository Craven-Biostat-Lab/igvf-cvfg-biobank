"""Helper constants and functions to maintain consistent entity notation."""

import pandas as pd

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