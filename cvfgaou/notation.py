"""Helper functions to maintain consistent entity notation."""


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